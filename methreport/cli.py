"""MethReport CLI — entry point for clinical methylation reporting."""

from __future__ import annotations

import logging
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
    help=(
        "Clinical-grade DNA methylation imprinting disorder reporter.\n\n"
        "Direct BAM input (phased or unphased) — recommended:\n\n"
        "  methreport run --bam sample.bam --ref t2t\n\n"
        "  Uses modkit internally for accurate strand-merging and HP splitting.\n"
        "  Falls back to internal MM/ML parser if modkit is not installed.\n\n"
        "Pre-computed BEDMethyl input (from modkit pileup):\n\n"
        "  methreport run --bed sample.unphased.bed --bed-hp1 sample.hp1.bed "
        "--bed-hp2 sample.hp2.bed --ref t2t"
    ),
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
    # Input — mutually exclusive: either BED or BAM
    bed: Optional[Path] = typer.Option(
        None, "--bed",
        help="[Recommended] modkit BEDMethyl file (unphased). Avoids MM/ML strand artefacts.",
    ),
    bed_hp1: Optional[Path] = typer.Option(
        None, "--bed-hp1",
        help="modkit BEDMethyl for haplotype 1 (optional, used with --bed).",
    ),
    bed_hp2: Optional[Path] = typer.Option(
        None, "--bed-hp2",
        help="modkit BEDMethyl for haplotype 2 (optional, used with --bed).",
    ),
    bam: Optional[Path] = typer.Option(
        None, "--bam",
        help="Indexed modbam file. Supports phased BAM (HP tags). Uses modkit internally if available.",
    ),
    call_threshold: float = typer.Option(
        0.5, "--call-threshold",
        help="(modbam only) MM/ML probability threshold for a methylation call.",
    ),
    # Shared options
    ref: str = typer.Option("t2t", "--ref", "-r", help="Reference genome: 't2t' or 'hg38'"),
    out: Path = typer.Option(Path("methreport_output"), "--out", "-o", help="Output directory"),
    sample_id: Optional[str] = typer.Option(None, "--sample", "-s", help="Sample ID (default: input filename stem)"),
    controls: Optional[Path] = typer.Option(None, "--controls", "-c", help="User-supplied controls TSV/XLSX"),
    replace_controls: bool = typer.Option(False, "--replace-controls", help="Replace (not supplement) bundled controls"),
    min_cov: int = typer.Option(5, "--min-cov", help="Minimum CpG coverage to include a site"),
    no_bed_tracks: bool = typer.Option(False, "--no-bed-tracks", help="Skip BED track output"),
    no_tsv: bool = typer.Option(False, "--no-tsv", help="Skip TSV export"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Debug logging"),
) -> None:
    """Run methylation imprinting analysis and generate HTML report."""
    _setup_logging(verbose)

    # --- Validate input mode ---
    if bed is None and bam is None:
        console.print("[red]Error:[/red] Provide either [cyan]--bed[/cyan] (recommended) or [cyan]--bam[/cyan].")
        raise typer.Exit(1)
    if bed is not None and bam is not None:
        console.print("[red]Error:[/red] [cyan]--bed[/cyan] and [cyan]--bam[/cyan] are mutually exclusive.")
        raise typer.Exit(1)

    use_bed = bed is not None

    if use_bed:
        if not bed.exists():
            console.print(f"[red]Error:[/red] BED file not found: {bed}")
            raise typer.Exit(1)
        if bed_hp1 and not bed_hp1.exists():
            console.print(f"[red]Error:[/red] HP1 BED not found: {bed_hp1}")
            raise typer.Exit(1)
        if bed_hp2 and not bed_hp2.exists():
            console.print(f"[red]Error:[/red] HP2 BED not found: {bed_hp2}")
            raise typer.Exit(1)
        input_display = str(bed)
    else:
        if not bam.exists():
            console.print(f"[red]Error:[/red] BAM not found: {bam}")
            raise typer.Exit(1)
        index_ok = Path(str(bam) + ".bai").exists() or Path(str(bam) + ".csi").exists()
        if not index_ok:
            console.print(
                f"[red]Error:[/red] BAM index not found. Run: [cyan]samtools index {bam}[/cyan]"
            )
            raise typer.Exit(1)
        input_display = str(bam)

    ref_norm = ref.lower()
    if ref_norm not in ("t2t", "hg38"):
        console.print(f"[red]Error:[/red] Unknown genome '{ref}'. Choose 't2t' or 'hg38'.")
        raise typer.Exit(1)

    input_path = bed if use_bed else bam
    sid = sample_id or input_path.stem
    out.mkdir(parents=True, exist_ok=True)

    from methreport import __version__
    mode_label = "BEDMethyl (modkit)" if use_bed else "modbam (MM/ML)"
    console.print(Panel(
        f"[bold blue]MethReport v{__version__}[/bold blue]\n"
        f"Sample: [cyan]{sid}[/cyan]  ·  Genome: [cyan]{ref_norm.upper()}[/cyan]  "
        f"·  Mode: [cyan]{mode_label}[/cyan]\n"
        f"Output: [cyan]{out}[/cyan]",
        title="MethReport",
        border_style="blue",
    ))

    if not use_bed:
        from methreport.modbam_pileup import modkit_version, MODKIT_MIN_VERSION
        mk_ver = modkit_version()
        if mk_ver and mk_ver >= MODKIT_MIN_VERSION:
            mk_str = ".".join(map(str, mk_ver))
            console.print(
                f"[dim]Using modkit v{mk_str} for strand-aware methylation extraction.[/dim]"
            )
        else:
            console.print(
                "[yellow]Note:[/yellow] modkit not found — using internal MM/ML parser with "
                "strand-normalisation fix. Install modkit for best accuracy."
            )

    # --- BAM validation (modbam mode only) ---
    bam_info: dict = {}
    if not use_bed:
        from methreport.reader import validate_bam
        with Progress(SpinnerColumn(), TextColumn("{task.description}"), console=console) as prog:
            t = prog.add_task("Validating BAM...", total=None)
            bam_info = validate_bam(bam)
            prog.remove_task(t)
        if not bam_info["has_mod_tags"]:
            console.print(
                "[yellow]Warning:[/yellow] No MM/ML modification tags detected in the first 1000 reads. "
                "Ensure this BAM was produced with basecaller modification calling enabled."
            )

    # --- Analysis ---
    from methreport.analysis import run_analysis

    with Progress(SpinnerColumn(), TextColumn("{task.description}"), console=console) as prog:
        t = prog.add_task("Analysing 14 imprinting DMRs...", total=None)
        analysis = run_analysis(
            genome=ref_norm,
            sample_id=sid,
            bed_unphased=bed,
            bed_hp1=bed_hp1,
            bed_hp2=bed_hp2,
            bam_path=bam,
            call_threshold=call_threshold,
            min_coverage=min_cov,
            user_controls_path=controls,
            replace_controls=replace_controls,
        )
        prog.remove_task(t)

    _print_summary(analysis)

    # --- Exports ---
    if not no_tsv:
        with Progress(SpinnerColumn(), TextColumn("{task.description}"), console=console) as prog:
            t = prog.add_task("Writing TSV exports...", total=None)
            from methreport.export import write_cpg_tsv, write_summary_tsv
            write_summary_tsv(analysis, out)
            write_cpg_tsv(analysis, out)
            prog.remove_task(t)

    if not no_bed_tracks:
        with Progress(SpinnerColumn(), TextColumn("{task.description}"), console=console) as prog:
            t = prog.add_task("Writing BED tracks...", total=None)
            from methreport.export import write_bed_track
            write_bed_track(analysis, out)
            prog.remove_task(t)

    with Progress(SpinnerColumn(), TextColumn("{task.description}"), console=console) as prog:
        t = prog.add_task("Generating HTML report...", total=None)
        from methreport.report import generate_report
        run_meta = {
            "Input": input_display,
            "Input mode": mode_label,
            "Min. coverage": str(min_cov),
            "Flagging method": "z-score vs controls" if True else "fixed threshold",
            "Controls": str(controls) if controls else "Bundled defaults",
        }
        if not use_bed:
            run_meta["Call threshold"] = str(call_threshold)
        report_path = generate_report(analysis, out, run_metadata=run_meta)
        prog.remove_task(t)

    console.print(f"\n[bold green]Done![/bold green]  Report → [cyan]{report_path}[/cyan]")

    if analysis.flagged_regions:
        console.print(f"\n[bold red]⚠ {len(analysis.flagged_regions)} region(s) flagged:[/bold red]")
        for r in analysis.flagged_regions:
            z_str = f"  z={r.z_score:+.2f}" if not (r.z_score != r.z_score) else ""
            console.print(
                f"  • {r.region.label}: {r.flag} ({r.mean_methylation:.1f}%){z_str}"
                f"  [{r.n_informative} informative CpGs]"
            )

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
    from methreport.modbam_pileup import bam_has_hp_tags, modkit_version, MODKIT_MIN_VERSION

    info = validate_bam(bam)
    has_hp = bam_has_hp_tags(bam)
    mk_ver = modkit_version()

    console.print(f"[bold]BAM:[/bold]          {info['path']}")
    console.print(f"[bold]Contigs:[/bold]      {info['n_contigs']} total, first: {', '.join(info['contigs'])}")
    console.print(f"[bold]Mod tags:[/bold]     {'[green]Yes — MM/ML present[/green]' if info['has_mod_tags'] else '[red]No MM/ML tags found[/red]'}")
    console.print(f"[bold]HP (phased):[/bold]  {'[green]Yes — HP tags detected[/green]' if has_hp else '[dim]No HP tags (unphased)[/dim]'}")

    if mk_ver and mk_ver >= MODKIT_MIN_VERSION:
        mk_str = ".".join(map(str, mk_ver))
        console.print(f"[bold]modkit:[/bold]       [green]v{mk_str} — will be used for strand-aware extraction[/green]")
    elif mk_ver:
        mk_str = ".".join(map(str, mk_ver))
        console.print(f"[bold]modkit:[/bold]       [yellow]v{mk_str} (< 0.2.0) — upgrade for best accuracy[/yellow]")
    else:
        console.print("[bold]modkit:[/bold]       [yellow]Not found — internal MM/ML parser will be used[/yellow]")


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
        table.add_row(r.name, r.label, r.disorder, r.chrom,
                      f"{r.start:,}", f"{r.end:,}", f"{r.end - r.start:,}")
    console.print(table)


def _print_summary(analysis) -> None:
    import numpy as np

    table = Table(title=f"Results — {analysis.sample_id}", border_style="blue", show_lines=True)
    table.add_column("Region", style="bold")
    table.add_column("Disorder")
    table.add_column("CpGs", justify="right")
    table.add_column("Info.", justify="right")
    table.add_column("Coverage", justify="right")
    table.add_column("Methylation", justify="right")
    table.add_column("HP1", justify="right")
    table.add_column("HP2", justify="right")
    table.add_column("z-score", justify="right")
    table.add_column("Status")

    for r in analysis.results:
        meth = r.mean_methylation
        hp1 = r.mean_methylation_hp1
        hp2 = r.mean_methylation_hp2
        z = r.z_score
        meth_str = f"{meth:.1f}%" if not np.isnan(meth) else "N/A"
        hp1_str  = f"{hp1:.1f}%"  if not np.isnan(hp1)  else "—"
        hp2_str  = f"{hp2:.1f}%"  if not np.isnan(hp2)  else "—"
        z_str    = f"{z:+.2f}"    if not np.isnan(z)     else "—"
        flag_style = "red bold" if r.flag in ("LOW", "HIGH") else "green"
        table.add_row(
            r.region.label, r.region.disorder,
            str(r.unphased.n_cpg), str(r.n_informative),
            f"{r.unphased.mean_coverage:.0f}×",
            meth_str, hp1_str, hp2_str, z_str,
            f"[{flag_style}]{r.flag}[/{flag_style}]",
        )

    console.print(table)


def main() -> None:
    app()


if __name__ == "__main__":
    main()
