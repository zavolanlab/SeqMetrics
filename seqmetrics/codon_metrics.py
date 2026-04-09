from collections import Counter
from typing import Dict, List
import math


def _relative_synonymous_codon_usage(
    codon_counts: Counter,
    syn_codons_by_aa: Dict[str, List[str]],
) -> Dict[str, float]:
    """Return RSCU per codon."""
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


def codon_adaptation_index_tissue(
    seq: str,
    tissue: str,
    tissue_codon_usage: Dict[str, Counter],
    syn_codons_by_aa: Dict[str, List[str]],
    min_weight: float = 0.01,
) -> float:
    """
    Tissue-specific CAI.

    Args:
        seq: CDS sequence (string).
        tissue: tissue name key in `tissue_codon_usage`.
        tissue_codon_usage: {tissue: Counter(codon -> count)}.
        syn_codons_by_aa: mapping aa -> list of synonymous codons.
        min_weight: floor for unknown codons / zero RSCU.

    Returns:
        CAI value for this CDS in the given tissue.
    """
    if tissue not in tissue_codon_usage:
        raise KeyError(f"Tissue '{tissue}' not found in tissue_codon_usage")

    ref_counts = tissue_codon_usage[tissue]
    rscu_ref = _relative_synonymous_codon_usage(ref_counts, syn_codons_by_aa)

    # relative adaptiveness weights per codon
    w: Dict[str, float] = {}
    for aa, codons in syn_codons_by_aa.items():
        rvals = [rscu_ref.get(c, 0.0) for c in codons]
        if not rvals:
            continue
        max_r = max(rvals)
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


