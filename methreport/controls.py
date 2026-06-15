"""Control data loading and reference range computation.

Bundled controls are derived from ctrls_2.xlsx / ctrls_hg38.xlsx in the
NanoImprint repository (carolinehey/NanoImprint).

Control data format (TSV, one row per CpG site per control sample):
    region_name  position  sample_id  methylation_pct

Users can supply additional or replacement controls via --controls flag.
"""

from __future__ import annotations

import logging
from importlib import resources
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

# Methylation range considered "imprinted" in controls (40-60% = heterozygous)
IMPRINT_MIN = 40.0
IMPRINT_MAX = 60.0

# Abnormal thresholds for fixed-threshold flagging (last-resort fallback only —
# z-score and phased allelic logic take priority when available).
# These are intentionally conservative to avoid false positives when controls
# and phasing are both absent. Coverage artefacts can lower the unphased mean
# to ~35% in normal samples, so we only flag at clearly extreme values.
FLAG_LOW = 20.0
FLAG_HIGH = 80.0


def load_bundled_controls(genome: str) -> pd.DataFrame:
    """Load the bundled control dataset for the given genome.

    Returns a DataFrame with columns:
        region_name, position, sample_id, methylation_pct
    """
    genome = genome.lower()
    filename = f"controls_{genome}.tsv"
    try:
        data_path = resources.files("methreport.data").joinpath(filename)
        df = pd.read_csv(str(data_path), sep="\t")
        log.debug("Loaded %d bundled control rows for %s", len(df), genome)
        return df
    except (FileNotFoundError, ModuleNotFoundError, TypeError):
        log.warning("Bundled controls for '%s' not found at %s — returning empty DataFrame", genome, filename)
        return pd.DataFrame(columns=["region_name", "position", "sample_id", "methylation_pct"])


def load_user_controls(path: Path | str) -> pd.DataFrame:
    """Load user-supplied controls TSV.

    Expected columns: region_name, position, sample_id, methylation_pct
    Missing columns cause a ValueError.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Controls file not found: {path}")

    suffix = path.suffix.lower()
    if suffix in (".xlsx", ".xls"):
        df = pd.read_excel(path)
    elif suffix in (".tsv", ".txt"):
        df = pd.read_csv(path, sep="\t")
    elif suffix == ".csv":
        df = pd.read_csv(path)
    else:
        raise ValueError(f"Unsupported controls format: {suffix}")

    required = {"region_name", "position", "sample_id", "methylation_pct"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Controls file missing columns: {missing}")

    log.info("Loaded %d user control rows from %s", len(df), path)
    return df[list(required)]


def merge_controls(
    genome: str,
    user_controls_path: Path | str | None = None,
    replace: bool = False,
) -> pd.DataFrame:
    """Return merged control DataFrame.

    If replace=True, bundled controls are replaced by user controls.
    Otherwise user controls are appended.
    """
    if replace and user_controls_path is not None:
        return load_user_controls(user_controls_path)

    bundled = load_bundled_controls(genome)
    if user_controls_path is not None:
        user = load_user_controls(user_controls_path)
        return pd.concat([bundled, user], ignore_index=True)
    return bundled


def compute_reference_ranges(
    controls: pd.DataFrame,
    region_name: str,
) -> pd.DataFrame:
    """Compute per-position mean ± SD from control data for a region.

    Returns DataFrame with columns:
        position, mean, sd, n_controls, lower_1sd, upper_1sd, lower_2sd, upper_2sd
    Only positions where controls show 40-60% methylation (imprinted pattern) are included.
    """
    subset = controls[controls["region_name"] == region_name].copy()
    if subset.empty:
        return pd.DataFrame(columns=["position", "mean", "sd", "n_controls",
                                     "lower_1sd", "upper_1sd", "lower_2sd", "upper_2sd"])

    # Identify imprinted positions: those where controls are 40-60%
    pos_stats = (
        subset.groupby("position")["methylation_pct"]
        .agg(mean="mean", sd="std", n_controls="count")
        .reset_index()
    )
    pos_stats["sd"] = pos_stats["sd"].fillna(0.0)
    pos_stats = pos_stats[
        pos_stats["mean"].between(IMPRINT_MIN, IMPRINT_MAX)
    ].copy()

    pos_stats["lower_1sd"] = pos_stats["mean"] - pos_stats["sd"]
    pos_stats["upper_1sd"] = pos_stats["mean"] + pos_stats["sd"]
    pos_stats["lower_2sd"] = pos_stats["mean"] - 2 * pos_stats["sd"]
    pos_stats["upper_2sd"] = pos_stats["mean"] + 2 * pos_stats["sd"]

    return pos_stats.sort_values("position").reset_index(drop=True)


def flag_result(mean_methylation: float) -> str:
    """Return a flag string using fixed thresholds (fallback when controls are absent).

    Prefer z_score_flag() when per-position control stats are available.
    """
    if np.isnan(mean_methylation):
        return "NA"
    if mean_methylation < FLAG_LOW:
        return "LOW"
    if mean_methylation > FLAG_HIGH:
        return "HIGH"
    return "NORMAL"


# Minimum number of informative CpGs required for a z-score flag to be trusted
MIN_INFORMATIVE_CPG = 3

# z-score magnitude that triggers a flag (~95th percentile under normal distribution)
Z_FLAG_THRESHOLD = 2.0

# Minimum control SD to avoid division by near-zero; corresponds to ~2% methylation
# measurement noise, which is realistic for nanopore at ≥10× coverage
_MIN_SD = 2.0


def z_score_flag(
    sample_df: pd.DataFrame,
    reference_df: pd.DataFrame,
    min_informative: int = MIN_INFORMATIVE_CPG,
) -> tuple[str, float, int]:
    """Flag a region using per-CpG z-scores against the control distribution.

    Only positions present in both the sample and the reference (i.e., positions
    where controls show an imprinted 40-60% pattern) are used. This avoids
    including uninformative CpGs that dilute the signal.

    Parameters
    ----------
    sample_df:
        DataFrame with columns [position, methylation_pct] from RegionMethylation.
    reference_df:
        Output of compute_reference_ranges() — columns [position, mean, sd, ...].
    min_informative:
        Minimum number of overlapping CpG positions needed to trust the z-score.
        Falls back to fixed-threshold flagging if fewer positions overlap.

    Returns
    -------
    (flag, mean_z, n_informative)
        flag          : "NORMAL" | "LOW" | "HIGH" | "NA"
        mean_z        : mean z-score across informative positions (nan if not computed)
        n_informative : number of positions used for the z-score calculation
    """
    if reference_df.empty or sample_df.empty:
        return "NA", float("nan"), 0

    merged = sample_df[["position", "methylation_pct"]].merge(
        reference_df[["position", "mean", "sd"]],
        on="position",
        how="inner",
    )
    n_informative = len(merged)

    if n_informative == 0:
        # No positional overlap at all — cannot make a meaningful comparison
        return "NA", float("nan"), 0

    if n_informative < min_informative:
        # Too few positions for a reliable z-score — use fixed threshold on what we have
        mean_meth = merged["methylation_pct"].mean()
        return flag_result(mean_meth), float("nan"), n_informative

    # Clip SD to a minimum to prevent extreme z-scores from near-zero control variance
    sd_clipped = merged["sd"].clip(lower=_MIN_SD)
    merged = merged.copy()
    merged["z"] = (merged["methylation_pct"] - merged["mean"]) / sd_clipped

    mean_z = float(merged["z"].mean())

    if mean_z < -Z_FLAG_THRESHOLD:
        return "LOW", mean_z, n_informative
    if mean_z > Z_FLAG_THRESHOLD:
        return "HIGH", mean_z, n_informative
    return "NORMAL", mean_z, n_informative
