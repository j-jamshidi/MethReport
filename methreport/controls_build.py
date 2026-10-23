"""Build an in-house control cohort from per-sample MethReport output.

Each analysed control sample leaves a `<sample>_cpg_methylation.tsv` behind. This
module pools those into the long control schema that `controls.py` accepts,

    region_name  position  sample_id  methylation_pct

keeping the unphased rows only, so the cohort yields a genuine per-position
mean, SD and control count. A provenance record is written beside it, because
a reference range is only as defensible as the account of how it was made.
"""

from __future__ import annotations

import json
import logging
import re
from datetime import date
from pathlib import Path

import pandas as pd

from methreport import __version__
from methreport.regions import get_regions

log = logging.getLogger(__name__)

_SUFFIX = "_cpg_methylation.tsv"
_REQUIRED = {"region", "haplotype", "position", "n_total", "methylation_pct"}

# Samples are identified in this repository by run number, never by episode
# accession. An id shaped like an accession is refused rather than recorded.
_ACCESSION = re.compile(r"^\d{2}[A-Z]\d{5,}")


def sample_id_from_path(path: Path) -> str:
    name = path.name
    return name[: -len(_SUFFIX)] if name.endswith(_SUFFIX) else path.stem


def build_controls(
    tsv_paths: list[Path],
    genome: str,
    min_cov: int = 5,
    sample_ids: list[str] | None = None,
) -> pd.DataFrame:
    """Pool per-CpG TSVs into the long control schema (unphased rows only)."""
    if not tsv_paths:
        raise ValueError("no per-CpG TSV files given")
    if sample_ids is not None and len(sample_ids) != len(tsv_paths):
        raise ValueError("sample_ids must match tsv_paths in length")

    label_to_name = {r.label: r.name for r in get_regions(genome)}
    ids = sample_ids or [sample_id_from_path(p) for p in tsv_paths]
    if len(set(ids)) != len(ids):
        raise ValueError("sample ids are not unique; pass them explicitly")
    for sid in ids:
        if _ACCESSION.match(sid):
            raise ValueError(
                f"sample id '{sid}' looks like an episode accession; use the "
                "laboratory run identifier instead"
            )

    frames = []
    for path, sid in zip(tsv_paths, ids):
        df = pd.read_csv(path, sep="\t")
        missing = _REQUIRED - set(df.columns)
        if missing:
            raise ValueError(f"{path.name}: missing columns {sorted(missing)}")

        df = df[(df["haplotype"] == "unphased") & (df["n_total"] >= min_cov)].copy()
        unknown = set(df["region"]) - set(label_to_name)
        if unknown:
            raise ValueError(
                f"{path.name}: regions {sorted(unknown)} are not defined for "
                f"'{genome}'; was this file made with a different --ref?"
            )
        df["region_name"] = df["region"].map(label_to_name)
        df["sample_id"] = sid
        frames.append(df[["region_name", "position", "sample_id", "methylation_pct"]])
        log.info("%s: %d unphased CpGs at >=%dx", sid, len(df), min_cov)

    out = pd.concat(frames, ignore_index=True)
    return out.sort_values(["region_name", "position", "sample_id"]).reset_index(drop=True)


def provenance(
    controls: pd.DataFrame,
    genome: str,
    min_cov: int,
    mod_code: str,
    basecaller_model: str,
    pipeline_version: str,
    notes: str = "",
) -> dict:
    """Describe the cohort: who is in it, how it was made, and n per position."""
    per_position = controls.groupby(["region_name", "position"])["sample_id"].nunique()
    per_region = {
        region: {
            "positions": int(len(grp)),
            "n_min": int(grp.min()),
            "n_median": float(grp.median()),
            "n_max": int(grp.max()),
        }
        for region, grp in per_position.groupby(level=0)
    }
    return {
        "created": date.today().isoformat(),
        "methreport_version": __version__,
        "genome": genome,
        "min_cov": min_cov,
        "modification": mod_code,
        "basecaller_model": basecaller_model,
        "upstream_pipeline_version": pipeline_version,
        "n_samples": int(controls["sample_id"].nunique()),
        "samples": sorted(controls["sample_id"].unique()),
        "n_positions": int(len(per_position)),
        "per_region": per_region,
        "notes": notes,
    }


def write_controls(controls: pd.DataFrame, record: dict, out: Path) -> tuple[Path, Path]:
    out.parent.mkdir(parents=True, exist_ok=True)
    controls.to_csv(out, sep="\t", index=False)
    prov = out.with_suffix(".provenance.json")
    prov.write_text(json.dumps(record, indent=2) + "\n")
    return out, prov
