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
        "Single sample:\n\n"
        "  methreport run --bam sample.bam --ref hg38\n\n"
        "Batch (one BAM per line in a text file):\n\n"
        "  methreport run --list samples.txt --ref hg38\n\n"
        "Pre-computed BEDMethyl input:\n\n"
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


def _parse_bam_list(path: Path) -> list[Path]:
    """Read a text file and return one BAM Path per non-blank, non-comment line."""
    bams = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        bams.append(Path(line))
    return bams


def _check_bam(bam: Path) -> str | None:
    """Return an error string if the BAM is missing or unindexed, else None."""
    if not bam.exists():
        return f"BAM not found: {bam}"
    if not (Path(str(bam) + ".bai").exists() or Path(str(bam) + ".csi").exists()):
        return f"BAM index not found (run: samtools index {bam})"
    return None


def _run_single_bam(
    bam: Path,
    ref_norm: str,
    out_dir: Path,
    sid: str,
    controls: Optional[Path],
    replace_controls: bool,
    min_cov: int,
    call_threshold: float,
    no_bed_tracks: bool,
    no_tsv: bool,
) -> tuple:
    """Run analysis on one BAM, write outputs to out_dir, return (analysis, report_path)."""
    from methreport.analysis import run_analysis
    from methreport.reader import validate_bam

    out_dir.mkdir(parents=True, exist_ok=True)

    bam_info = validate_bam(bam)
    if not bam_info["has_mod_tags"]:
        console.print(
            f"[yellow]Warning:[/yellow] No MM/ML tags detected in {bam.name}. "
            "Ensure modification calling was enabled during basecalling."
        )

    analysis = run_analysis(
        genome=ref_norm,
        sample_id=sid,
        bam_path=bam,
        call_threshold=call_threshold,
        min_coverage=min_cov,
        user_controls_path=controls,
        replace_controls=replace_controls,
    )

    if not no_tsv:
        from methreport.export import write_cpg_tsv, write_summary_tsv
        write_summary_tsv(analysis, out_dir)
        write_cpg_tsv(analysis, out_dir)

    if not no_bed_tracks:
        from methreport.export import write_bed_track
        write_bed_track(analysis, out_dir)

    from methreport.report import generate_report
    run_meta = {
        "Input": str(bam),
        "Input mode": "modbam (MM/ML)",
        "Min. coverage": str(min_cov),
        "Call threshold": str(call_threshold),
        "Controls": str(controls) if controls else "Bundled defaults",
    }
    report_path = generate_report(analysis, out_dir, run_metadata=run_meta)
    return analysis, report_path


@app.command()
def run(
    # Input — mutually exclusive groups
    bam: Optional[Path] = typer.Option(
        None, "--bam",
        help="Indexed modbam file (single sample). Supports phased BAM (HP tags).",
    ),
    bam_list: Optional[Path] = typer.Option(
        None, "--list",
        help="Text file with one BAM path per line. Runs batch analysis on all. "
             "Each sample gets its own subdirectory under --out.",
    ),
    bed: Optional[Path] = typer.Option(
        None, "--bed",
        help="modkit BEDMethyl file (unphased). Single-sample alternative to --bam.",
    ),
    bed_hp1: Optional[Path] = typer.Option(
        None, "--bed-hp1",
        help="modkit BEDMethyl for haplotype 1 (used with --bed).",
    ),
    bed_hp2: Optional[Path] = typer.Option(
        None, "--bed-hp2",
        help="modkit BEDMethyl for haplotype 2 (used with --bed).",
    ),
    call_threshold: float = typer.Option(
        0.5, "--call-threshold",
        help="MM/ML probability threshold for methylation calls (modbam only).",
    ),
    # Shared options
    ref: str = typer.Option("t2t", "--ref", "-r", help="Reference genome: 't2t' or 'hg38'"),
    out: Path = typer.Option(Path("methreport_output"), "--out", "-o", help="Output directory"),
    sample_id: Optional[str] = typer.Option(
        None, "--sample", "-s",
        help="Sample ID (single-sample mode only; ignored with --list).",
    ),
    controls: Optional[Path] = typer.Option(None, "--controls", "-c", help="User-supplied controls TSV/XLSX"),
    replace_controls: bool = typer.Option(False, "--replace-controls", help="Replace (not supplement) bundled controls"),
    min_cov: int = typer.Option(5, "--min-cov", help="Minimum CpG coverage to include a site"),
    no_bed_tracks: bool = typer.Option(False, "--no-bed-tracks", help="Skip BED track output"),
    no_tsv: bool = typer.Option(False, "--no-tsv", help="Skip TSV export"),
    verbose: bool = typer.Option(False, "--verbose", "-v", help="Debug logging"),
) -> None:
    """Run methylation imprinting analysis and generate HTML report(s)."""
    _setup_logging(verbose)

    # --- Validate input mode ---
    n_inputs = sum(x is not None for x in (bam, bam_list, bed))
    if n_inputs == 0:
        console.print(
            "[red]Error:[/red] Provide one of [cyan]--bam[/cyan], [cyan]--list[/cyan], or [cyan]--bed[/cyan]."
        )
        raise typer.Exit(1)
    if n_inputs > 1:
        console.print("[red]Error:[/red] [cyan]--bam[/cyan], [cyan]--list[/cyan], and [cyan]--bed[/cyan] are mutually exclusive.")
        raise typer.Exit(1)

    ref_norm = ref.lower()
    if ref_norm not in ("t2t", "hg38"):
        console.print(f"[red]Error:[/red] Unknown genome '{ref}'. Choose 't2t' or 'hg38'.")
        raise typer.Exit(1)

    from methreport import __version__
    from methreport.modbam_pileup import MODKIT_MIN_VERSION, modkit_version

    # ══════════════════════════════════════════════════════════════════════════
    # BATCH MODE  (--list)
    # ══════════════════════════════════════════════════════════════════════════
    if bam_list is not None:
        if not bam_list.exists():
            console.print(f"[red]Error:[/red] List file not found: {bam_list}")
            raise typer.Exit(1)

        bam_paths = _parse_bam_list(bam_list)
        if not bam_paths:
            console.print(f"[red]Error:[/red] No BAM paths found in {bam_list}")
            raise typer.Exit(1)

        mk_ver = modkit_version()
        mk_label = (
            f"modkit v{'.'.join(map(str, mk_ver))}"
            if mk_ver and mk_ver >= MODKIT_MIN_VERSION
            else "internal MM/ML parser"
        )
        out.mkdir(parents=True, exist_ok=True)

        console.print(Panel(
            f"[bold blue]MethReport v{__version__} — Batch Mode[/bold blue]\n"
            f"Samples: [cyan]{len(bam_paths)}[/cyan]  ·  Genome: [cyan]{ref_norm.upper()}[/cyan]  "
            f"·  Extractor: [cyan]{mk_label}[/cyan]\n"
            f"Output root: [cyan]{out}[/cyan]",
            title="MethReport Batch",
            border_style="blue",
        ))

        # Pre-flight checks
        errors = []
        for bp in bam_paths:
            err = _check_bam(bp)
            if err:
                errors.append(err)
        if errors:
            console.print(f"[red]Pre-flight errors ({len(errors)}):[/red]")
            for e in errors:
                console.print(f"  [red]✗[/red] {e}")
            console.print("Fix the above issues and re-run.")
            raise typer.Exit(1)

        all_analyses = []
        failed = []

        for i, bp in enumerate(bam_paths, 1):
            sid = bp.stem
            sample_out = out / sid
            console.rule(f"[bold]{i}/{len(bam_paths)}: {sid}[/bold]")
            try:
                with Progress(SpinnerColumn(), TextColumn("{task.description}"), console=console) as prog:
                    t = prog.add_task(f"Analysing {sid}...", total=None)
                    analysis, report_path = _run_single_bam(
                        bam=bp,
                        ref_norm=ref_norm,
                        out_dir=sample_out,
                        sid=sid,
                        controls=controls,
                        replace_controls=replace_controls,
                        min_cov=min_cov,
                        call_threshold=call_threshold,
                        no_bed_tracks=no_bed_tracks,
                        no_tsv=no_tsv,
                    )
                    prog.remove_task(t)

                flagged = analysis.flagged_regions
                if flagged:
                    flag_str = ", ".join(f"{r.region.label} ({r.flag})" for r in flagged)
                    console.print(f"[red]  ⚠ {len(flagged)} flagged:[/red] {flag_str}")
                else:
                    console.print("[green]  ✓ All regions normal[/green]")
                console.print(f"  Report → [cyan]{report_path}[/cyan]")
                all_analyses.append(analysis)

            except Exception as exc:
                console.print(f"[red]  ✗ FAILED: {exc}[/red]")
                failed.append((sid, str(exc)))

        # Batch summary
        _print_batch_summary(all_analyses)
        _write_batch_summary_tsv(all_analyses, out)

        if failed:
            console.print(f"\n[red bold]{len(failed)} sample(s) failed:[/red bold]")
            for sid, err in failed:
                console.print(f"  [red]✗[/red] {sid}: {err}")
            raise typer.Exit(1)

        console.print(
            f"\n[bold green]Batch complete![/bold green]  "
            f"{len(all_analyses)}/{len(bam_paths)} samples processed.  "
            f"Summary → [cyan]{out / 'batch_summary.tsv'}[/cyan]"
        )
        raise typer.Exit(0)

    # ══════════════════════════════════════════════════════════════════════════
    # SINGLE SAMPLE — BAM mode
    # ══════════════════════════════════════════════════════════════════════════
    if bam is not None:
        err = _check_bam(bam)
        if err:
            console.print(f"[red]Error:[/red] {err}")
            raise typer.Exit(1)

        sid = sample_id or bam.stem
        out.mkdir(parents=True, exist_ok=True)

        mk_ver = modkit_version()
        mk_label = (
            f"modkit v{'.'.join(map(str, mk_ver))}"
            if mk_ver and mk_ver >= MODKIT_MIN_VERSION
            else "internal MM/ML parser"
        )

        console.print(Panel(
            f"[bold blue]MethReport v{__version__}[/bold blue]\n"
            f"Sample: [cyan]{sid}[/cyan]  ·  Genome: [cyan]{ref_norm.upper()}[/cyan]  "
            f"·  Extractor: [cyan]{mk_label}[/cyan]\n"
            f"Output: [cyan]{out}[/cyan]",
            title="MethReport",
            border_style="blue",
        ))

        with Progress(SpinnerColumn(), TextColumn("{task.description}"), console=console) as prog:
            t = prog.add_task("Analysing 14 imprinting DMRs...", total=None)
            analysis, report_path = _run_single_bam(
                bam=bam,
                ref_norm=ref_norm,
                out_dir=out,
                sid=sid,
                controls=controls,
                replace_controls=replace_controls,
                min_cov=min_cov,
                call_threshold=call_threshold,
                no_bed_tracks=no_bed_tracks,
                no_tsv=no_tsv,
            )
            prog.remove_task(t)

        _print_summary(analysis)
        console.print(f"\n[bold green]Done![/bold green]  Report → [cyan]{report_path}[/cyan]")

        if analysis.flagged_regions:
            console.print(f"\n[bold red]⚠ {len(analysis.flagged_regions)} region(s) flagged:[/bold red]")
            for r in analysis.flagged_regions:
                z_str = f"  z={r.z_score:+.2f}" if r.z_score == r.z_score else ""
                console.print(
                    f"  • {r.region.label}: {r.flag} ({r.mean_methylation:.1f}%){z_str}"
                    f"  [{r.n_informative} informative CpGs]"
                )
        raise typer.Exit(0)

    # ══════════════════════════════════════════════════════════════════════════
    # SINGLE SAMPLE — BED mode
    # ══════════════════════════════════════════════════════════════════════════
    if bed is not None:
        if not bed.exists():
            console.print(f"[red]Error:[/red] BED file not found: {bed}")
            raise typer.Exit(1)
        if bed_hp1 and not bed_hp1.exists():
            console.print(f"[red]Error:[/red] HP1 BED not found: {bed_hp1}")
            raise typer.Exit(1)
        if bed_hp2 and not bed_hp2.exists():
            console.print(f"[red]Error:[/red] HP2 BED not found: {bed_hp2}")
            raise typer.Exit(1)

        sid = sample_id or bed.stem
        out.mkdir(parents=True, exist_ok=True)

        console.print(Panel(
            f"[bold blue]MethReport v{__version__}[/bold blue]\n"
            f"Sample: [cyan]{sid}[/cyan]  ·  Genome: [cyan]{ref_norm.upper()}[/cyan]  "
            f"·  Mode: [cyan]BEDMethyl[/cyan]\n"
            f"Output: [cyan]{out}[/cyan]",
            title="MethReport",
            border_style="blue",
        ))

        from methreport.analysis import run_analysis

        with Progress(SpinnerColumn(), TextColumn("{task.description}"), console=console) as prog:
            t = prog.add_task("Analysing 14 imprinting DMRs...", total=None)
            analysis = run_analysis(
                genome=ref_norm,
                sample_id=sid,
                bed_unphased=bed,
                bed_hp1=bed_hp1,
                bed_hp2=bed_hp2,
                min_coverage=min_cov,
                user_controls_path=controls,
                replace_controls=replace_controls,
            )
            prog.remove_task(t)

        _print_summary(analysis)

        if not no_tsv:
            from methreport.export import write_cpg_tsv, write_summary_tsv
            write_summary_tsv(analysis, out)
            write_cpg_tsv(analysis, out)

        if not no_bed_tracks:
            from methreport.export import write_bed_track
            write_bed_track(analysis, out)

        from methreport.report import generate_report
        run_meta = {
            "Input": str(bed),
            "Input mode": "BEDMethyl (modkit)",
            "Min. coverage": str(min_cov),
            "Controls": str(controls) if controls else "Bundled defaults",
        }
        report_path = generate_report(analysis, out, run_metadata=run_meta)
        console.print(f"\n[bold green]Done![/bold green]  Report → [cyan]{report_path}[/cyan]")

        if analysis.flagged_regions:
            console.print(f"\n[bold red]⚠ {len(analysis.flagged_regions)} region(s) flagged:[/bold red]")
            for r in analysis.flagged_regions:
                z_str = f"  z={r.z_score:+.2f}" if r.z_score == r.z_score else ""
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

    from methreport.modbam_pileup import MODKIT_MIN_VERSION, bam_has_hp_tags, modkit_version
    from methreport.reader import validate_bam

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
    table.add_column("Note")
    for r in regions:
        note = "[yellow]⚠ unreliable[/yellow]" if r.unreliable else ""
        table.add_row(r.name, r.label, r.disorder, r.chrom,
                      f"{r.start:,}", f"{r.end:,}", f"{r.end - r.start:,}", note)
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
        if r.flag in ("LOW", "HIGH"):
            flag_cell = f"[red bold]{r.flag}[/red bold]"
        elif r.flag == "UNRELIABLE":
            flag_cell = "[yellow]UNRELIABLE ⚠[/yellow]"
        else:
            flag_cell = f"[green]{r.flag}[/green]"
        table.add_row(
            r.region.label, r.region.disorder,
            str(r.unphased.n_cpg), str(r.n_informative),
            f"{r.unphased.mean_coverage:.0f}×",
            meth_str, hp1_str, hp2_str, z_str,
            flag_cell,
        )

    console.print(table)


def _print_batch_summary(analyses: list) -> None:
    """Print a condensed one-row-per-sample table after batch run."""
    import numpy as np

    console.rule("[bold]Batch Summary[/bold]")
    table = Table(border_style="blue", show_lines=True)
    table.add_column("Sample", style="bold")
    table.add_column("Flagged", justify="center")
    table.add_column("Flagged regions")

    for analysis in analyses:
        flagged = analysis.flagged_regions
        n = len(flagged)
        if n:
            flag_detail = ", ".join(
                f"{r.region.label} ({r.flag})" for r in flagged
            )
            status_str = f"[red bold]{n}[/red bold]"
        else:
            flag_detail = "[dim]—[/dim]"
            status_str = "[green]0[/green]"
        table.add_row(analysis.sample_id, status_str, flag_detail)

    console.print(table)


def _write_batch_summary_tsv(analyses: list, out_dir: Path) -> None:
    """Write a combined TSV with one row per region per sample."""
    import pandas as pd

    rows = []
    for analysis in analyses:
        df = analysis.summary_table()
        df.insert(0, "sample_id", analysis.sample_id)
        rows.append(df)

    if not rows:
        return

    combined = pd.concat(rows, ignore_index=True)
    tsv_path = out_dir / "batch_summary.tsv"
    combined.to_csv(tsv_path, sep="\t", index=False)
    import logging
    logging.getLogger(__name__).info("Batch summary TSV → %s", tsv_path)


def main() -> None:
    app()


if __name__ == "__main__":
    main()
