"""modbam reader — extracts per-CpG 5mC methylation from BAM files with MM/ML tags.

Supports:
  - Unphased reads (all reads combined)
  - Phased reads split by HP:i:1 / HP:i:2 haplotag
  - Direct modbam input (Oxford Nanopore modkit-compatible BAM)

The MM/ML tag spec (SAM v1.7 §1.7) encodes modified bases as a string of
delta-encoded positions in the MM tag and probabilities in the ML tag.
We interpret ML values as methylation probability: prob = ML / 255.
A read-level call is made at threshold 0.5 (configurable).
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd
import pysam

from methreport.regions import DMRegion

log = logging.getLogger(__name__)

Haplotype = Literal["unphased", "hp1", "hp2"]

# Minimum ML probability to call a base as methylated
DEFAULT_CALL_THRESHOLD = 0.5
# Minimum reads covering a CpG to include it in analysis
DEFAULT_MIN_COVERAGE = 5


@dataclass
class CpGSite:
    chrom: str
    position: int       # 0-based genomic position of the C in CpG
    n_total: int = 0
    n_methyl: int = 0

    @property
    def methylation_pct(self) -> float:
        if self.n_total == 0:
            return float("nan")
        return 100.0 * self.n_methyl / self.n_total


@dataclass
class RegionMethylation:
    region: DMRegion
    haplotype: Haplotype
    sites: list[CpGSite] = field(default_factory=list)

    @property
    def mean_coverage(self) -> float:
        if not self.sites:
            return 0.0
        return float(np.mean([s.n_total for s in self.sites]))

    @property
    def mean_methylation_pct(self) -> float:
        vals = [s.methylation_pct for s in self.sites if not np.isnan(s.methylation_pct)]
        if not vals:
            return float("nan")
        return float(np.mean(vals))

    @property
    def n_cpg(self) -> int:
        return len(self.sites)

    def to_dataframe(self) -> pd.DataFrame:
        if not self.sites:
            return pd.DataFrame(columns=["chrom", "position", "n_total", "n_methyl", "methylation_pct"])
        rows = [
            {
                "chrom": s.chrom,
                "position": s.position,
                "n_total": s.n_total,
                "n_methyl": s.n_methyl,
                "methylation_pct": s.methylation_pct,
            }
            for s in self.sites
        ]
        return pd.DataFrame(rows)


def _parse_mm_ml_tags(
    read: pysam.AlignedSegment,
    call_threshold: float,
) -> dict[int, bool]:
    """Return {ref_pos: is_methylated} for all 5mC calls on a read.

    Handles both forward and reverse strand CpGs.
    Returns an empty dict if MM/ML tags are absent or malformed.
    """
    if not read.has_tag("MM") or not read.has_tag("ML"):
        return {}

    mm_str: str = read.get_tag("MM")
    ml_array: list[int] = list(read.get_tag("ML"))

    # Parse MM tag — we only want C+m (5mC) entries
    # MM format: <base><strand><mod>[,<delta>]*;...
    seq = read.query_sequence
    if seq is None:
        return {}
    seq_bytes = seq.encode()

    # Build mapping from query position → reference position
    # pysam gives us (query_pos, ref_pos) pairs via get_aligned_pairs
    q2r: dict[int, int] = {}
    for qp, rp in read.get_aligned_pairs(matches_only=True):
        if rp is not None:
            q2r[qp] = rp

    result: dict[int, bool] = {}
    ml_offset = 0

    for mod_entry in mm_str.rstrip(";").split(";"):
        if not mod_entry:
            continue
        parts = mod_entry.split(",")
        header = parts[0]  # e.g. "C+m" or "C+m?"

        # Only process 5mC on cytosines
        if not (header.startswith("C+m") or header.startswith("C-m")):
            # Skip non-5mC modifications but advance ML offset
            # Count how many positions this mod entry covers
            n_positions = len(parts) - 1
            ml_offset += n_positions
            continue

        strand_char = header[1]  # '+' forward, '-' reverse
        deltas = [int(x) for x in parts[1:] if x]
        n_positions = len(deltas)

        # Walk through the read sequence finding cytosine positions
        if strand_char == "+":
            search_base = ord("C")
        else:
            search_base = ord("G")

        cyto_positions = [i for i, b in enumerate(seq_bytes) if b == search_base]

        cursor = -1
        for i, delta in enumerate(deltas):
            cursor += delta + 1
            if cursor >= len(cyto_positions):
                break
            q_pos = cyto_positions[cursor]
            r_pos = q2r.get(q_pos)
            if r_pos is None:
                continue
            prob = ml_array[ml_offset + i] / 255.0
            result[r_pos] = prob >= call_threshold

        ml_offset += n_positions

    return result


def extract_region_methylation(
    bam_path: Path | str,
    region: DMRegion,
    call_threshold: float = DEFAULT_CALL_THRESHOLD,
    min_coverage: int = DEFAULT_MIN_COVERAGE,
) -> dict[Haplotype, RegionMethylation]:
    """Extract per-CpG methylation for a DMR from a modbam file.

    Returns a dict with keys "unphased", "hp1", "hp2".
    "hp1" and "hp2" are empty if the BAM is not phased.
    """
    bam_path = Path(bam_path)
    if not bam_path.exists():
        raise FileNotFoundError(f"BAM file not found: {bam_path}")

    index_path = Path(str(bam_path) + ".bai")
    if not index_path.exists():
        csi_path = Path(str(bam_path) + ".csi")
        if not csi_path.exists():
            raise FileNotFoundError(
                f"BAM index not found. Run: samtools index {bam_path}"
            )

    # Accumulate counts: {haplotype: {position: [n_total, n_methyl]}}
    counts: dict[Haplotype, dict[int, list[int]]] = {
        "unphased": {},
        "hp1": {},
        "hp2": {},
    }

    chrom = region.chrom
    # Remove 'chr' prefix if contig names in BAM don't use it
    bam = pysam.AlignmentFile(str(bam_path), "rb")

    # Detect whether the BAM uses 'chr' prefixes
    contigs = {c["SN"] for c in bam.header.to_dict().get("SQ", [])}
    if chrom not in contigs:
        alt = chrom.replace("chr", "") if chrom.startswith("chr") else f"chr{chrom}"
        if alt in contigs:
            chrom = alt
        else:
            log.warning("Contig %s not found in BAM header, skipping region %s", chrom, region.name)
            bam.close()
            return _empty_result(region)

    try:
        fetch_iter = bam.fetch(chrom, region.start, region.end)
    except ValueError as exc:
        log.warning("Could not fetch region %s: %s", region.name, exc)
        bam.close()
        return _empty_result(region)

    reads_processed = 0
    reads_with_mods = 0

    for read in fetch_iter:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        reads_processed += 1

        hp: Haplotype = "unphased"
        if read.has_tag("HP"):
            hp_val = int(read.get_tag("HP"))
            if hp_val == 1:
                hp = "hp1"
            elif hp_val == 2:
                hp = "hp2"

        mod_calls = _parse_mm_ml_tags(read, call_threshold)
        if not mod_calls:
            continue
        reads_with_mods += 1

        for ref_pos, is_methyl in mod_calls.items():
            if not (region.start <= ref_pos < region.end):
                continue
            for h in ("unphased", hp) if hp != "unphased" else ("unphased",):
                if ref_pos not in counts[h]:
                    counts[h][ref_pos] = [0, 0]
                counts[h][ref_pos][0] += 1
                counts[h][ref_pos][1] += int(is_methyl)

    bam.close()

    if reads_processed > 0:
        log.debug(
            "Region %s: %d reads processed, %d with mod tags (%.0f%%)",
            region.name,
            reads_processed,
            reads_with_mods,
            100 * reads_with_mods / reads_processed,
        )

    result: dict[Haplotype, RegionMethylation] = {}
    for hp in ("unphased", "hp1", "hp2"):
        sites = []
        for pos in sorted(counts[hp]):
            n_total, n_methyl = counts[hp][pos]
            if n_total >= min_coverage:
                sites.append(CpGSite(chrom=region.chrom, position=pos, n_total=n_total, n_methyl=n_methyl))
        result[hp] = RegionMethylation(region=region, haplotype=hp, sites=sites)

    return result


def _empty_result(region: DMRegion) -> dict[Haplotype, RegionMethylation]:
    return {
        "unphased": RegionMethylation(region=region, haplotype="unphased"),
        "hp1": RegionMethylation(region=region, haplotype="hp1"),
        "hp2": RegionMethylation(region=region, haplotype="hp2"),
    }


def validate_bam(bam_path: Path | str) -> dict:
    """Return basic QC info about the BAM file."""
    bam_path = Path(bam_path)
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    header = bam.header.to_dict()
    contigs = [sq["SN"] for sq in header.get("SQ", [])]
    has_mod_reads = False
    n_checked = 0
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped:
            continue
        if read.has_tag("MM"):
            has_mod_reads = True
            break
        n_checked += 1
        if n_checked > 1000:
            break
    bam.close()
    return {
        "path": str(bam_path),
        "contigs": contigs[:10],
        "n_contigs": len(contigs),
        "has_mod_tags": has_mod_reads,
    }
