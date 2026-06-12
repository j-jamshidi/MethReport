"""MethReport CLI — entry point for clinical methylation reporting."""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.logging import RichHandler
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich.table import Table

app = typer.Typer(
    name="methreport",
    help="Clinical-grade DNA methylation imprinting disorder reporter.",
    add_completion=False,
    pretty_exceptions_enable=False,
)

console = Console()


def _setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(message)s",
        handlers=[RichHandler(console=console, rich_tracebacks=True, show_path=False)],
    )


@app.command()
def run(
    bam: Path = typer.Argument(..., help="Indexed modbam file (.bam)"),
    ref: str = typer.Option("t2t", "--ref", "-r", help="Reference genome: 't2t' or 'hg38'"),
    out: Path = typer.Option(Path("methreport_output"), "--out", "-o", help="Output directory"),
    sample_id: Optional[str] = typer.Option(None, "--sample", "-s", help="Sample identifier (default: BAM filename)"),
    controls: Optional[Path] = typer.Option(None, "--controls", "-c", help="User-supplied controls TSV/XLSX"),
    replace_controls: bool = typer.Option(False, "--replace-controls", help="Replace (not supplement) bundled controls"),
    min_cov: int = typer.Option(5, "--min-cov", help="Minimum CpG coverage to include a site"),
    call_threshold: float = typer.Option(0.5, "--call-threshold", help="MM/ML probability threshold for methylation call"),
    no_bed: bool = typer.Option(False, "--no-bed", help="Skip BED track output"),
    no_tsv: bool = typer.Option(False, "--no-tsv", help="Skip TSV export"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Debug logging"),
) -> None:
    """Analyse a modbam file and generate a methylation imprinting report."""
    _setup_logging(verbose)
    log = logging.getLogger(__name__)

    # --- Validate inputs ---
    if not bam.exists():
        console.print(f"[red]Error:[/red] BAM file not found: {bam}")
        raise typer.Exit(1)

    index_ok = Path(str(bam) + ".bai").exists() or Path(str(bam) + ".csi").exists()
    if not index_ok:
        console.print(f"[red]Error:[/red] BAM index not found. Run: [cyan]samtools index {bam}[/cyan]")
        raise typer.Exit(1)

    ref_norm = ref.lower()
    if ref_norm not in ("t2t", "hg38"):
        console.print(f"[red]Error:[/red] Unknown genome '{ref}'. Choose 't2t' or 'hg38'.")
        raise typer.Exit(1)

    sid = sample_id or bam.stem
    out.mkdir(parents=True, exist_ok=True)

    # --- Header ---
    from methreport import __version__
    console.print(Panel(
        f"[bold blue]MethReport v{__version__}[/bold blue]\n"
        f"Sample: [cyan]{sid}[/cyan]  |  Genome: [cyan]{ref_norm.upper()}[/cyan]  |  Output: [cyan]{out}[/cyan]",
        title="MethReport",
        border_style="blue",
    ))

    # --- BAM validation ---
    from methreport.reader import validate_bam
    with Progress(SpinnerColumn(), TextColumn("{task.description}"), console=console) as prog:
        t = prog.add_task("Validating BAM...", total=None)
        bam_info = validate_bam(bam)
        prog.remove_task(t)

    if not bam_info["has_mod_tags"]:
        console.print("[yellow]Warning:[/yellow] No MM/ML modification tags found in the first 1000 reads. "
                      "Ensure this is a modbam file produced with --mod enabled.")

    # --- Analysis ---
    from methreport.analysis import run_analysis

    with Progress(SpinnerColumn(), TextColumn("{task.description}"), console=console) as prog:
        t = prog.add_task(f"Analysing {len(_get_regions(ref_norm))} DMRs...", total=None)
        analysis = run_analysis(
            bam_path=bam,
            genome=ref_norm,
            sample_id=sid,
            user_controls_path=controls,
            replace_controls=replace_controls,
            call_threshold=call_threshold,
            min_coverage=min_cov,
        )
        prog.remove_task(t)

    # --- Print summary ---
    _print_summary(analysis)

    # --- Export ---
    from methreport.export import write_all_exports
    from methreport.report import generate_report

    exports = {}
    if not no_tsv:
        with Progress(SpinnerColumn(), TextColumn("{task.description}"), console=console) as prog:
            t = prog.add_task("Writing TSV exports...", total=None)
            from methreport.export import write_cpg_tsv, write_summary_tsv
            exports["summary_tsv"] = write_summary_tsv(analysis, out)
            exports["cpg_tsv"] = write_cpg_tsv(analysis, out)
            prog.remove_task(t)

    if not no_bed:
        with Progress(SpinnerColumn(), TextColumn("{task.description}"), console=console) as prog:
            t = prog.add_task("Writing BED tracks...", total=None)
            from methreport.export import write_bed_track
            exports["bed_tracks"] = write_bed_track(analysis, out)
            prog.remove_task(t)

    with Progress(SpinnerColumn(), TextColumn("{task.description}"), console=console) as prog:
        t = prog.add_task("Generating HTML report...", total=None)
        run_meta = {
            "BAM file": str(bam),
            "BAM contigs": str(bam_info["n_contigs"]),
            "Min. coverage": str(min_cov),
            "Call threshold": str(call_threshold),
            "Controls": str(controls) if controls else "Bundled defaults",
        }
        report_path = generate_report(analysis, out, run_metadata=run_meta)
        prog.remove_task(t)

    # --- Done ---
    console.print(f"\n[bold green]Done![/bold green]  Report: [cyan]{report_path}[/cyan]")

    if analysis.flagged_regions:
        console.print(f"[bold red]{len(analysis.flagged_regions)} region(s) flagged:[/bold red]")
        for r in analysis.flagged_regions:
            console.print(f"  • {r.region.label}: {r.flag} ({r.mean_methylation:.1f}%)")

    raise typer.Exit(0)


@app.command()
def validate(
    bam: Path = typer.Argument(..., help="BAM file to inspect"),
) -> None:
    """Quick check: verify a BAM has modbam MM/ML tags and print contig info."""
    _setup_logging(False)
    if not bam.exists():
        console.print(f"[red]Not found:[/red] {bam}")
        raise typer.Exit(1)

    from methreport.reader import validate_bam
    info = validate_bam(bam)
    console.print(f"[bold]BAM:[/bold]  {info['path']}")
    console.print(f"[bold]Contigs:[/bold] {info['n_contigs']}")
    console.print(f"[bold]First contigs:[/bold] {', '.join(info['contigs'])}")
    console.print(f"[bold]Has mod tags:[/bold] {'[green]Yes[/green]' if info['has_mod_tags'] else '[red]No[/red]'}")


@app.command("list-regions")
def list_regions(
    ref: str = typer.Option("t2t", "--ref", "-r", help="Reference genome: 't2t' or 'hg38'"),
) -> None:
    """List all DMR regions for the given reference genome."""
    _setup_logging(False)
    from methreport.regions import get_regions
    regions = get_regions(ref.lower())

    table = Table(title=f"DMR Regions ({ref.upper()})", border_style="blue")
    table.add_column("Name", style="cyan")
    table.add_column("Label")
    table.add_column("Disorder")
    table.add_column("Chrom")
    table.add_column("Start", justify="right")
    table.add_column("End", justify="right")
    table.add_column("Size (bp)", justify="right")

    for r in regions:
        table.add_row(
            r.name,
            r.label,
            r.disorder,
            r.chrom,
            f"{r.start:,}",
            f"{r.end:,}",
            f"{r.end - r.start:,}",
        )

    console.print(table)


def _get_regions(genome: str) -> list:
    from methreport.regions import get_regions
    return get_regions(genome)


def _print_summary(analysis) -> None:
    from methreport.controls import FLAG_HIGH, FLAG_LOW

    table = Table(title=f"Results — {analysis.sample_id}", border_style="blue", show_lines=True)
    table.add_column("Region", style="bold")
    table.add_column("Disorder")
    table.add_column("CpGs", justify="right")
    table.add_column("Coverage", justify="right")
    table.add_column("Methylation", justify="right")
    table.add_column("HP1", justify="right")
    table.add_column("HP2", justify="right")
    table.add_column("Status")

    import numpy as np

    for r in analysis.results:
        meth = r.mean_methylation
        meth_str = f"{meth:.1f}%" if not np.isnan(meth) else "N/A"
        hp1 = r.mean_methylation_hp1
        hp2 = r.mean_methylation_hp2
        hp1_str = f"{hp1:.1f}%" if not np.isnan(hp1) else "—"
        hp2_str = f"{hp2:.1f}%" if not np.isnan(hp2) else "—"
        status_style = "red bold" if r.flag in ("LOW", "HIGH") else "green"
        flag_str = f"[{status_style}]{r.flag}[/{status_style}]"
        table.add_row(
            r.region.label,
            r.region.disorder,
            str(r.unphased.n_cpg),
            f"{r.unphased.mean_coverage:.0f}×",
            meth_str,
            hp1_str,
            hp2_str,
            flag_str,
        )

    console.print(table)


def main() -> None:
    app()


if __name__ == "__main__":
    main()
