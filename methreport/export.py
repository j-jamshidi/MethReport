"""Export analysis results to TSV tables and IGV-ready BED tracks."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from methreport.analysis import SampleAnalysis
from methreport.reader import RegionMethylation

log = logging.getLogger(__name__)


def write_summary_tsv(analysis: SampleAnalysis, out_dir: Path) -> Path:
    """Write per-region summary table as TSV."""
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{analysis.sample_id}_summary.tsv"
    analysis.summary_table().to_csv(path, sep="\t", index=False)
    log.info("Summary TSV → %s", path)
    return path


def write_cpg_tsv(analysis: SampleAnalysis, out_dir: Path) -> Path:
    """Write per-CpG methylation values for all regions (unphased + phased)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{analysis.sample_id}_cpg_methylation.tsv"

    frames = []
    for result in analysis.results:
        for hp_key in ("unphased", "hp1", "hp2"):
            rm: RegionMethylation = getattr(result, hp_key) if hp_key == "unphased" else getattr(result, hp_key)
            if hp_key == "unphased":
                rm = result.unphased
            elif hp_key == "hp1":
                rm = result.hp1
            else:
                rm = result.hp2

            if rm.n_cpg == 0:
                continue
            df = rm.to_dataframe()
            df.insert(0, "region", result.region.label)
            df.insert(1, "haplotype", hp_key)
            frames.append(df)

    if frames:
        pd.concat(frames, ignore_index=True).to_csv(path, sep="\t", index=False)
    else:
        pd.DataFrame(columns=["region", "haplotype", "chrom", "position",
                               "n_total", "n_methyl", "methylation_pct"]).to_csv(path, sep="\t", index=False)

    log.info("Per-CpG TSV → %s", path)
    return path


def write_bed_track(analysis: SampleAnalysis, out_dir: Path) -> list[Path]:
    """Write BED9 tracks suitable for IGV / UCSC browser.

    Produces three files:
      - <sample>_unphased.bed
      - <sample>_hp1.bed   (empty if not phased)
      - <sample>_hp2.bed   (empty if not phased)

    Score column (0-1000) encodes methylation percentage.
    itemRgb encodes methylation level: blue=low, white=mid, red=high.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []

    for hp_key in ("unphased", "hp1", "hp2"):
        lines = [
            f'track name="{analysis.sample_id} ({hp_key})" '
            f'description="MethReport 5mC methylation" '
            f'itemRgb="On" useScore=0\n'
        ]

        for result in analysis.results:
            if hp_key == "unphased":
                rm = result.unphased
            elif hp_key == "hp1":
                rm = result.hp1
            else:
                rm = result.hp2

            for site in rm.sites:
                pct = site.methylation_pct
                score = min(1000, int(pct * 10))
                r, g, b = _methylation_rgb(pct)
                lines.append(
                    f"{site.chrom}\t{site.position}\t{site.position + 1}\t"
                    f"CpG\t{score}\t.\t{site.position}\t{site.position + 1}\t"
                    f"{r},{g},{b}\n"
                )

        out_path = out_dir / f"{analysis.sample_id}_{hp_key}.bed"
        out_path.write_text("".join(lines))
        log.info("BED track → %s", out_path)
        paths.append(out_path)

    return paths


def _methylation_rgb(pct: float) -> tuple[int, int, int]:
    """Map methylation % to an RGB color: blue (0%) → white (50%) → red (100%)."""
    t = max(0.0, min(1.0, pct / 100.0))
    if t <= 0.5:
        s = t / 0.5
        r = int(255 * s)
        g = int(255 * s)
        b = 255
    else:
        s = (t - 0.5) / 0.5
        r = 255
        g = int(255 * (1.0 - s))
        b = int(255 * (1.0 - s))
    return r, g, b


def write_all_exports(analysis: SampleAnalysis, out_dir: Path) -> dict[str, Path | list[Path]]:
    """Write all export files and return paths dict."""
    return {
        "summary_tsv": write_summary_tsv(analysis, out_dir),
        "cpg_tsv": write_cpg_tsv(analysis, out_dir),
        "bed_tracks": write_bed_track(analysis, out_dir),
    }
