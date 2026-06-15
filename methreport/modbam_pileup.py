"""Extract methylation from phased modbam files without pre-generating BED files.

Strategy (tried in order):
  1. modkit subprocess — if modkit ≥ 0.2.0 is on PATH.
     Runs two commands over all 14 regions at once:
       a) unphased pileup  → combines + and − strand CpG calls correctly
       b) partitioned pileup (--partition-tag HP) → auto-splits HP1 / HP2
  2. Internal pysam fallback — if modkit is absent.
     Applies the strand-normalisation fix directly: reverse-strand C calls are
     shifted by −1 to align with the + strand CpG coordinate (same logic as
     modkit --combine-strands).

Why modkit is preferred
-----------------------
Our MM/ML parser finds C positions in each read's query sequence and maps them
to reference coordinates via get_aligned_pairs(). For a CpG at ref position p:

  • Forward-strand read  →  C at ref pos p   (correct)
  • Reverse-strand read  →  C at ref pos p+1 (complement base — off by one)

Without normalisation the same CpG generates counts at both p and p+1, halving
apparent coverage and creating adjacent noise spikes. modkit handles this with
--combine-strands; our fallback applies the −1 shift manually.

Reference
---------
modkit documentation: https://nanoporetech.github.io/modkit/
MM/ML spec: SAMv1 §1.7
"""

from __future__ import annotations

import logging
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
import pysam

from methreport.reader import (
    CpGSite,
    RegionMethylation,
    _empty_result,
    _parse_mm_ml_tags,
)
from methreport.regions import DMRegion, get_regions

log = logging.getLogger(__name__)

MODKIT_MIN_VERSION = (0, 2, 0)

# ─────────────────────────────────────────────────────────────────────────────
# modkit helpers
# ─────────────────────────────────────────────────────────────────────────────

def modkit_path() -> str | None:
    return shutil.which("modkit")


def modkit_version() -> tuple[int, ...] | None:
    """Parse `modkit --version` output → (major, minor, patch) or None."""
    path = modkit_path()
    if not path:
        return None
    try:
        r = subprocess.run([path, "--version"], capture_output=True, text=True, timeout=10)
        for line in r.stdout.splitlines():
            parts = line.strip().split()
            if len(parts) >= 2 and parts[0].lower() == "modkit":
                return tuple(int(x) for x in parts[1].split(".")[:3])
    except Exception:
        pass
    return None


def modkit_available() -> bool:
    v = modkit_version()
    return v is not None and v >= MODKIT_MIN_VERSION


# ─────────────────────────────────────────────────────────────────────────────
# BAM introspection
# ─────────────────────────────────────────────────────────────────────────────

def bam_has_hp_tags(bam_path: Path, scan_limit: int = 2000) -> bool:
    """Return True if any mapped read in the BAM carries an HP tag."""
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    n = 0
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped:
            continue
        if read.has_tag("HP"):
            bam.close()
            return True
        n += 1
        if n >= scan_limit:
            break
    bam.close()
    return False


def _bam_uses_chr_prefix(bam_path: Path) -> bool:
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    contigs = [sq["SN"] for sq in bam.header.to_dict().get("SQ", [])]
    bam.close()
    return any(c.startswith("chr") for c in contigs[:10])


# ─────────────────────────────────────────────────────────────────────────────
# Region BED writer (for modkit --include-bed)
# ─────────────────────────────────────────────────────────────────────────────

def _write_region_bed(genome: str, out_path: Path, use_chr: bool) -> None:
    """Write a temporary BED file covering all 14 DMR windows."""
    with open(out_path, "w") as fh:
        for r in get_regions(genome):
            chrom = r.chrom if use_chr else r.chrom.lstrip("chr")
            fh.write(f"{chrom}\t{r.start}\t{r.end}\t{r.name}\n")


# ─────────────────────────────────────────────────────────────────────────────
# modkit subprocess runner
# ─────────────────────────────────────────────────────────────────────────────

def _run_modkit(cmd: list[str], label: str, timeout: int = 300) -> bool:
    log.debug("modkit %s: %s", label, " ".join(cmd))
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        if r.returncode != 0:
            log.warning("modkit %s exited %d\nstderr: %s", label, r.returncode, r.stderr[-1500:])
            return False
        return True
    except subprocess.TimeoutExpired:
        log.warning("modkit %s timed out after %ds", label, timeout)
        return False
    except Exception as exc:
        log.warning("modkit %s error: %s", label, exc)
        return False


def _find_partition_file(search_dir: Path, hp_value: str) -> Path | None:
    """Locate the modkit partition output file for the given HP value.

    modkit names partition outputs as {prefix}_{value}.bed.
    We search flexibly to handle minor naming differences across versions.
    """
    # Exact match first
    for p in search_dir.glob(f"*_{hp_value}.bed"):
        if "ungrouped" not in p.name:
            return p
    # Looser match
    for p in search_dir.glob("*.bed"):
        if "ungrouped" not in p.name and hp_value in p.stem.split("_"):
            return p
    return None


# ─────────────────────────────────────────────────────────────────────────────
# modkit extraction
# ─────────────────────────────────────────────────────────────────────────────

def _extract_via_modkit(
    bam_path: Path,
    genome: str,
    min_coverage: int,
    call_threshold: float,
) -> dict[str, dict[str, RegionMethylation]]:
    """Run modkit pileup once for all regions, return per-region results."""
    from methreport.bed_reader import _load_bed, extract_region_from_bed

    mk = modkit_path()
    regions = get_regions(genome)
    use_chr = _bam_uses_chr_prefix(bam_path)
    has_hp = bam_has_hp_tags(bam_path)

    with tempfile.TemporaryDirectory(prefix="methreport_") as tmp:
        tmpdir = Path(tmp)
        region_bed = tmpdir / "regions.bed"
        _write_region_bed(genome, region_bed, use_chr=use_chr)

        # Shared flags
        # Note: modkit 0.6+ has no --no-header flag; default output has no header.
        # Our _load_bed already skips any '#'-prefixed header lines defensively.
        base = [
            mk, "pileup", str(bam_path),
            "--filter-threshold", str(call_threshold),
            "--combine-strands",   # merge + / − strand CpG calls
            "--cpg",               # restrict to CpG dinucleotides
            "--include-bed", str(region_bed),
        ]

        # ── Unphased ──
        unphased_bed = tmpdir / "unphased.bed"
        if not _run_modkit(base + [str(unphased_bed)], "unphased pileup"):
            raise RuntimeError("modkit unphased pileup failed")
        unphased_df = _load_bed(unphased_bed)
        log.info("modkit unphased: %d 5mC sites across all regions", len(unphased_df))

        # ── Phased ──
        empty_df = pd.DataFrame(columns=["chrom", "position", "n_valid_cov", "pct_modified"])
        hp1_df, hp2_df = empty_df.copy(), empty_df.copy()

        if has_hp:
            phase_prefix = str(tmpdir / "hp")
            ok = _run_modkit(
                base + [phase_prefix, "--partition-tag", "HP"],
                "phased pileup",
            )
            if ok:
                p1 = _find_partition_file(tmpdir, "1")
                p2 = _find_partition_file(tmpdir, "2")
                if p1:
                    hp1_df = _load_bed(p1)
                    log.info("modkit HP1: %d sites", len(hp1_df))
                if p2:
                    hp2_df = _load_bed(p2)
                    log.info("modkit HP2: %d sites", len(hp2_df))
                if not p1 and not p2:
                    log.warning(
                        "modkit partition output not found in %s — "
                        "BAM may not be haplotagged despite HP reads detected", tmpdir
                    )
        else:
            log.info("No HP tags detected — skipping phased pileup")

        # ── Build result dict ──
        results: dict[str, dict[str, RegionMethylation]] = {}
        for region in regions:
            results[region.name] = {
                "unphased": extract_region_from_bed(unphased_df, region, "unphased", min_coverage),
                "hp1":      extract_region_from_bed(hp1_df,     region, "hp1",      min_coverage),
                "hp2":      extract_region_from_bed(hp2_df,     region, "hp2",      min_coverage),
            }
    return results


# ─────────────────────────────────────────────────────────────────────────────
# pysam fallback with strand-normalisation fix
# ─────────────────────────────────────────────────────────────────────────────

def _extract_via_pysam(
    bam_path: Path,
    genome: str,
    min_coverage: int,
    call_threshold: float,
) -> dict[str, dict[str, RegionMethylation]]:
    """Strand-aware MM/ML parser — corrects reverse-strand CpG position by −1.

    For a CpG at + strand position p:
      forward-strand read → MM/ML reports ref pos p       (correct)
      reverse-strand read → MM/ML reports ref pos p+1    (off by one)

    Fix: subtract 1 from all ref positions on reverse-strand reads.
    This merges both strands to position p, matching modkit --combine-strands.
    """
    regions = get_regions(genome)
    results: dict[str, dict[str, RegionMethylation]] = {
        r.name: _empty_result(r) for r in regions
    }

    # counts[region_name][haplotype][position] = [n_total, n_methyl]
    Counts = dict[str, dict[str, dict[int, list[int]]]]
    counts: Counts = {
        r.name: {"unphased": {}, "hp1": {}, "hp2": {}}
        for r in regions
    }

    bam = pysam.AlignmentFile(str(bam_path), "rb")
    contigs = {c["SN"] for c in bam.header.to_dict().get("SQ", [])}

    for region in regions:
        chrom = region.chrom
        if chrom not in contigs:
            alt = chrom.lstrip("chr") if chrom.startswith("chr") else f"chr{chrom}"
            if alt in contigs:
                chrom = alt
            else:
                log.warning("Contig %s not in BAM, skipping %s", chrom, region.name)
                continue

        try:
            reads = bam.fetch(chrom, region.start, region.end)
        except ValueError as exc:
            log.warning("Could not fetch %s: %s", region.name, exc)
            continue

        n_reads = n_with_mods = 0
        for read in reads:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            n_reads += 1

            hp = "unphased"
            if read.has_tag("HP"):
                v = int(read.get_tag("HP"))
                hp = "hp1" if v == 1 else "hp2" if v == 2 else "unphased"

            mod_calls = _parse_mm_ml_tags(read, call_threshold)
            if not mod_calls:
                continue
            n_with_mods += 1

            is_rev = read.is_reverse
            for raw_pos, is_methyl in mod_calls.items():
                # Strand normalisation: reverse-strand C → shift to CpG + strand pos
                ref_pos = raw_pos - 1 if is_rev else raw_pos

                if not (region.start <= ref_pos < region.end):
                    continue

                targets = ["unphased"] if hp == "unphased" else ["unphased", hp]
                for h in targets:
                    c = counts[region.name][h]
                    if ref_pos not in c:
                        c[ref_pos] = [0, 0]
                    c[ref_pos][0] += 1
                    c[ref_pos][1] += int(is_methyl)

        if n_reads > 0:
            log.debug(
                "%s: %d reads, %d with mods (%.0f%%)",
                region.name, n_reads, n_with_mods, 100 * n_with_mods / n_reads,
            )

    bam.close()

    for region in regions:
        for hp in ("unphased", "hp1", "hp2"):
            sites = [
                CpGSite(chrom=region.chrom, position=pos,
                        n_total=c[0], n_methyl=c[1])
                for pos, c in sorted(counts[region.name][hp].items())
                if c[0] >= min_coverage
            ]
            results[region.name][hp] = RegionMethylation(
                region=region, haplotype=hp, sites=sites
            )

    return results


# ─────────────────────────────────────────────────────────────────────────────
# Public entry point
# ─────────────────────────────────────────────────────────────────────────────

def extract_from_bam(
    bam_path: Path | str,
    genome: str,
    min_coverage: int,
    call_threshold: float,
) -> dict[str, dict[str, RegionMethylation]]:
    """Extract per-region methylation from a (phased) modbam file.

    Automatically uses modkit if available, else falls back to the internal
    strand-aware pysam parser.

    Returns
    -------
    {region_name: {"unphased": RegionMethylation, "hp1": ..., "hp2": ...}}
    """
    bam_path = Path(bam_path)

    mk_ver = modkit_version()
    if mk_ver and mk_ver >= MODKIT_MIN_VERSION:
        log.info(
            "Using modkit %s for methylation extraction (--combine-strands --cpg)",
            ".".join(map(str, mk_ver)),
        )
        try:
            return _extract_via_modkit(bam_path, genome, min_coverage, call_threshold)
        except Exception as exc:
            log.warning("modkit extraction failed: %s — falling back to internal parser", exc)
    else:
        if mk_ver:
            log.warning(
                "modkit %s < 0.2.0, falling back to internal parser. "
                "Upgrade modkit for best accuracy.",
                ".".join(map(str, mk_ver)),
            )
        else:
            log.info(
                "modkit not found. Using internal MM/ML parser with strand-normalisation fix. "
                "Install modkit for best accuracy: https://github.com/nanoporetech/modkit"
            )

    return _extract_via_pysam(bam_path, genome, min_coverage, call_threshold)
