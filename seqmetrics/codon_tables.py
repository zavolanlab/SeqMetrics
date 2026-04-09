# seqmetrics/codon_tables.py

from __future__ import annotations
from collections import Counter
from typing import Dict, List, Set, Tuple
import math

########################################
# Basic nucleotide and codon utilities
########################################

def gc_at_content(seq: str) -> Tuple[float, float, float, float]:
    """
    Return nucleotide composition (A, G, C, T) as fractions of non-N bases.
    """
    seq = seq.upper()
    clean = seq.replace("N", "")
    length = len(clean) or 1
    a = clean.count("A") / length
    t = clean.count("T") / length
    g = clean.count("G") / length
    c = clean.count("C") / length
    return a, g, c, t

def codon_counts(seq: str) -> Counter:
    """
    Return codon composition for an in-frame CDS (no stop trimmed).
    """
    seq = seq.upper()
    codons = [seq[i:i+3] for i in range(0, len(seq), 3) if len(seq[i:i+3]) == 3]
    return Counter(codons)

########################################
# Tissue-specific codon usage
########################################

def load_tissue_codon_usage(path: str) -> Dict[str, Counter]:
    """
    Load tissue-specific codon usage from a table with a header.

    Expected format (tab- or space-separated):

        Tissue  AAA  AAC  AAG  ...  TTT
        brain   10   5    7    ...  3
        liver   20   7    9    ...  1
        ...

    - First column: tissue name (string).
    - Remaining columns: each column header is a codon (e.g. AAA, AAC,...).
    - Each subsequent line: tissue name + integer counts.

    Returns:
        dict: {tissue_name: Counter({codon: count, ...}), ...}
    """
    tissue_tables: Dict[str, Counter] = {}

    with open(path, "r", encoding="utf-8") as fh:
        header = fh.readline()
        if not header:
            raise ValueError("Codon-usage file is empty")

        # Split header into column names; first must be "Tissue" (or similar)
        header_cols = header.strip().split()
        if len(header_cols) < 2:
            raise ValueError("Header must contain at least 'Tissue' and one codon")

        tissue_col = header_cols[0]
        codon_cols = [c.upper() for c in header_cols[1:]]

        # basic sanity check: enforce 3-letter codons made of ACGT (optional)
        for c in codon_cols:
            if len(c) != 3 or any(b not in "ACGT" for b in c):
                raise ValueError(
                    f"Header column '{c}' does not look like a codon; "
                    "please check the input file format."
                )

        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            if len(parts) != len(header_cols):
                raise ValueError(
                    f"Line has {len(parts)} fields but header has {len(header_cols)}: {line}"
                )

            tissue = parts[0]
            counts_str = parts[1:]

            counter = Counter()
            for codon, val in zip(codon_cols, counts_str):
                # support both plain integers and things like "123," if CSV was converted badly
                val_clean = val.rstrip(",")
                counter[codon] += int(float(val_clean))

            tissue_tables[tissue] = counter

    return tissue_tables

########################################
# RSCU, optimal codons, and CAI per tissue
########################################

def _rscu_from_counts(
    codon_counts: Counter,
    syn_codons_by_aa: Dict[str, List[str]],
) -> Dict[str, float]:
    """
    Compute RSCU (Relative Synonymous Codon Usage) per codon
    from codon counts and synonymous groups.
    """
    rscu: Dict[str, float] = {}
    for aa, codons in syn_codons_by_aa.items():
        total = sum(codon_counts.get(c, 0) for c in codons)
        k = len(codons)
        if total == 0:
            for c in codons:
                rscu[c] = 0.0
        else:
            for c in codons:
                rscu[c] = (codon_counts.get(c, 0) * k) / total
    return rscu

def optimal_codons_for_tissue(
    tissue: str,
    tissue_codon_usage: Dict[str, Counter],
    syn_codons_by_aa: Dict[str, List[str]],
    rscu_threshold: float = 0.9,
) -> Set[str]:
    """
    Determine optimal codons for a tissue based on its codon usage.

    For each amino acid, codons whose RSCU is >= rscu_threshold * max_RSCU
    among its synonymous set are considered "optimal".
    """
    if tissue not in tissue_codon_usage:
        raise KeyError(f"Tissue '{tissue}' not found in tissue_codon_usage")

    counts = tissue_codon_usage[tissue]
    rscu = _rscu_from_counts(counts, syn_codons_by_aa)

    optimal: Set[str] = set()
    for aa, codons in syn_codons_by_aa.items():
        rscu_vals = [rscu.get(c, 0.0) for c in codons]
        if not rscu_vals:
            continue
        max_r = max(rscu_vals)
        if max_r == 0:
            continue
        for c in codons:
            if rscu.get(c, 0.0) >= rscu_threshold * max_r:
                optimal.add(c)

    return optimal

def codon_adaptation_index_tissue(
    seq: str,
    tissue: str,
    tissue_codon_usage: Dict[str, Counter],
    syn_codons_by_aa: Dict[str, List[str]],
    min_weight: float = 0.01,
) -> float:
    """
    Tissue-specific Codon Adaptation Index (CAI).

    Uses tissue-specific codon counts as the reference to compute RSCU
    and relative adaptiveness weights w_i. CAI is the geometric mean
    of w_i along the sequence.
    """
    if tissue not in tissue_codon_usage:
        raise KeyError(f"Tissue '{tissue}' not found in tissue_codon_usage")

    ref_counts = tissue_codon_usage[tissue]
    rscu_ref = _rscu_from_counts(ref_counts, syn_codons_by_aa)

    # Relative adaptiveness weights
    w: Dict[str, float] = {}
    for aa, codons in syn_codons_by_aa.items():
        rscu_vals = [rscu_ref.get(c, 0.0) for c in codons]
        if not rscu_vals:
            continue
        max_r = max(rscu_vals)
        if max_r == 0:
            for c in codons:
                w[c] = min_weight
            continue
        for c in codons:
            w[c] = rscu_ref.get(c, 0.0) / max_r

    codons_seq = [seq[i:i+3].upper() for i in range(0, len(seq), 3)
                  if len(seq[i:i+3]) == 3]
    if not codons_seq:
        return float("nan")

    logs = []
    for c in codons_seq:
        wi = w.get(c, min_weight)
        if wi <= 0:
            wi = min_weight
        logs.append(math.log(wi))

    return math.exp(sum(logs) / len(logs))

def fraction_optimal_codons_tissue(
    seq: str,
    tissue: str,
    tissue_optimal_codons: Dict[str, Set[str]],
) -> float:
    """
    Fraction of codons in the CDS that are optimal for a given tissue.
    """
    if tissue not in tissue_optimal_codons:
        raise KeyError(f"Tissue '{tissue}' not found in tissue_optimal_codons")

    optimal = tissue_optimal_codons[tissue]

    codons_seq = [seq[i:i+3].upper() for i in range(0, len(seq), 3)
                  if len(seq[i:i+3]) == 3]
    if not codons_seq:
        return 0.0

    opt_count = sum(1 for c in codons_seq if c in optimal)
    return opt_count / len(codons_seq)

