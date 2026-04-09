# trna_weights.py

from __future__ import annotations
from dataclasses import dataclass
from collections import Counter, defaultdict
from typing import Dict, List, Tuple, Optional
import math
import csv
import os

# ----------------------------------------------------------------------
# Standard nuclear genetic code (human) codon -> amino acid (1-letter)
# ----------------------------------------------------------------------

CODON_TO_AA: Dict[str, str] = {
    # Phe
    "TTT": "F", "TTC": "F",
    # Leu
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    # Ile
    "ATT": "I", "ATC": "I", "ATA": "I",
    # Met
    "ATG": "M",
    # Val
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    # Ser
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    # Pro
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    # Thr
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    # Ala
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    # Tyr
    "TAT": "Y", "TAC": "Y",
    # His
    "CAT": "H", "CAC": "H",
    # Gln
    "CAA": "Q", "CAG": "Q",
    # Asn
    "AAT": "N", "AAC": "N",
    # Lys
    "AAA": "K", "AAG": "K",
    # Asp
    "GAT": "D", "GAC": "D",
    # Glu
    "GAA": "E", "GAG": "E",
    # Cys
    "TGT": "C", "TGC": "C",
    # Trp
    "TGG": "W",
    # Arg
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    # Gly
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    # Stops
    "TAA": "*", "TAG": "*", "TGA": "*",
}

# ----------------------------------------------------------------------
# Wobble pairing scores (codon 3rd base, anticodon 1st base, in RNA)
# ----------------------------------------------------------------------

# codon_third_base : anticodon_base -> score
PAIRING_SCORES: Dict[Tuple[str, str], float] = {
    ("G", "U"): 0.56,
    ("U", "G"): 0.68,
    ("A", "U"): 1.0,
    ("U", "A"): 1.0,
    ("G", "C"): 1.0,
    ("C", "G"): 1.0,
    ("I", "C"): 0.99,
    ("I", "A"): 0.28,
    # inferred from “scores calculated based on the rules above”
    ("I", "U"): 1.0,
}

# Watson–Crick pairs (codon_base, anticodon_base) in RNA
WC_PAIRS = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")}


def _dna_to_rna(s: str) -> str:
    return s.upper().replace("T", "U")


@dataclass
class TRNAWeightTable:
    """
    Stores:
      - tGCN: amino acid -> anticodon -> gene copy number
      - W_by_codon: raw adaptation (sum_j s_ij * tGCN_ij)
      - w_by_codon: normalized weights W_i / W_max (or mean(w) if W_i=0)
    """
    species: str
    tGCN: Dict[str, Dict[str, int]]
    W_by_codon: Dict[str, float]
    w_by_codon: Dict[str, float]

    # ------------------------------------------------------------------
    # 1) Build from GtRNAdb-like tRNA gene file
    # ------------------------------------------------------------------
    @classmethod
    def from_trna_gene_file(
        cls,
        species: str,
        path: str,
        anticodon_col: str = "AntiCodon",
        aa_col: str = "Type",       # or "Isotype" if that matches better
        note_col: str = "Note",
        high_conf_phrase: str = "high confidence set",
    ) -> "TRNAWeightTable":
        """
        Parse a tRNA gene table (tab-delimited) like:

        SeqName tRNA# Begin End Type AntiCodon ... Note
        chr1    4     ...   ... Gly  CCC       ... high confidence set

        and tabulate tGCN[aa][anticodon] = count of high-confidence genes.
        Then build codon–anticodon interaction scores and codon weights.
        """
        if not os.path.exists(path):
            raise FileNotFoundError(path)

        # 1) parse tRNA gene table, build tGCN
        tGCN: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))

        with open(path, "r", encoding="utf-8") as fh:
            reader = csv.reader(fh, delimiter="\t")
            header = None
            for row in reader:
                if not row or all(not x.strip() for x in row):
                    continue
                if header is None:
                    header = [h.strip() for h in row]
                    continue
                # skip separator lines if present
                if set("".join(row).strip()) <= set("-"):
                    continue

                values = {h: (row[i].strip() if i < len(row) else "")
                          for i, h in enumerate(header)}

                note = values.get(note_col, "")
                if high_conf_phrase and high_conf_phrase not in note:
                    continue  # keep only high-confidence

                anticodon = values.get(anticodon_col, "")
                aa_name = values.get(aa_col, "")

                if not anticodon or not aa_name:
                    continue

                # anticodon as RNA 5'->3'
                ac_rna = _dna_to_rna(anticodon)

                # crude mapping Type/Isotype (3-letter) -> one-letter AA where possible
                aa_code = _aa_to_one_letter(aa_name)
                if aa_code is None:
                    continue

                tGCN[aa_code][ac_rna] += 1

        # 2) build W_i and w_i
        W_by_codon, w_by_codon = _build_codon_weights(tGCN)

        return cls(
            species=species,
            tGCN={aa: dict(ac_counts) for aa, ac_counts in tGCN.items()},
            W_by_codon=W_by_codon,
            w_by_codon=w_by_codon,
        )

    # ------------------------------------------------------------------
    # 2) Compute tAI for a sequence (exclude first/last codon)
    # ------------------------------------------------------------------
    def tAI(self, codons_in_seq: List[str]) -> float:
        """
        tAI = geometric mean of internal codon weights (excluding first and last codon).

        Uses precomputed self.w_by_codon.
        """
        if len(codons_in_seq) <= 2:
            return float("nan")

        internal_codons = [c.upper() for c in codons_in_seq[1:-1]]
        weights = []
        for c in internal_codons:
            w = self.w_by_codon.get(c)
            if w is None:
                # if codon not seen at all, approximate by mean weight
                if self.w_by_codon:
                    w = sum(self.w_by_codon.values()) / len(self.w_by_codon)
                else:
                    return float("nan")
            if w <= 0:
                continue
            weights.append(w)

        if not weights:
            return float("nan")

        log_sum = sum(math.log(w) for w in weights)
        return math.exp(log_sum / len(weights))


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def _aa_to_one_letter(name: str) -> Optional[str]:
    """
    Map GtRNAdb 'Type'/'Isotype' values to one-letter AA codes.
    E.g. 'Gly' -> 'G', 'Glu' -> 'E', 'Ile' -> 'I', etc.
    """
    name = name.strip()
    # common 3-letter forms
    MAP = {
        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D",
        "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
        "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
        "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
        "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
        "Sup": "*", "Sec": "U",
    }
    if name in MAP:
        return MAP[name]
    # if already single letter
    if len(name) == 1 and name.isalpha():
        return name.upper()
    return None


def _pairing_score_for(codon: str, anticodon: str) -> float:
    """
    Score s_ij for a given codon–anticodon pair (both 5'->3', DNA or RNA).
    Returns 0 if they cannot pair according to the rules.
    """
    codon_rna = _dna_to_rna(codon)
    ac_rna = _dna_to_rna(anticodon)

    if len(codon_rna) != 3 or len(ac_rna) != 3:
        return 0.0

    c0, c1, c2 = codon_rna[0], codon_rna[1], codon_rna[2]
    a0, a1, a2 = ac_rna[0], ac_rna[1], ac_rna[2]

    # positions 1 and 2 (codon[0]/anticodon[2], codon[1]/anticodon[1]) must be Watson-Crick
    if (c0, a2) not in WC_PAIRS or (c1, a1) not in WC_PAIRS:
        return 0.0

    # wobble position: codon[2] vs anticodon[0]
    # anticodons starting with A are effectively inosine (I) at position 34
    if a0 == "A":
        a0_eff = "I"
    else:
        a0_eff = a0

    key = (c2, a0_eff)
    return PAIRING_SCORES.get(key, 0.0)


def _build_codon_weights(
    tGCN: Dict[str, Dict[str, int]],
) -> Tuple[Dict[str, float], Dict[str, float]]:
    """
    From tGCN[aa][anticodon] and codon–AA mapping, build:

      W_i = sum_j s_ij * tGCN_ij
      w_i = W_i / max(W) if W_i != 0, else mean(w)

    Returns:
      (W_by_codon, w_by_codon)
    """
    # 1) W_i
    W_by_codon: Dict[str, float] = {}

    # Build list of all (aa, codon) we care about (exclude stops)
    aa_codons: Dict[str, List[str]] = defaultdict(list)
    for codon, aa in CODON_TO_AA.items():
        if aa == "*":
            continue
        aa_codons[aa].append(codon)

    for aa, codons in aa_codons.items():
        anticodons_for_aa = tGCN.get(aa, {})
        if not anticodons_for_aa:
            # if no tRNAs annotated for this AA, set W_i = 0 for all its codons
            for codon in codons:
                W_by_codon[codon] = 0.0
            continue

        for codon in codons:
            W_i = 0.0
            for anticodon, copies in anticodons_for_aa.items():
                s_ij = _pairing_score_for(codon, anticodon)
                if s_ij <= 0:
                    continue
                W_i += s_ij * copies
            W_by_codon[codon] = W_i

    # 2) Normalize to w_i
    W_values = [v for v in W_by_codon.values() if v > 0]
    if not W_values:
        # degenerate case: all zero
        return W_by_codon, {c: 0.0 for c in W_by_codon}

    W_max = max(W_values)
    w_by_codon: Dict[str, float] = {}
    non_zero_w: List[float] = []

    # first pass: codons with W_i > 0
    for codon, W_i in W_by_codon.items():
        if W_i > 0:
            w = W_i / W_max
            w_by_codon[codon] = w
            non_zero_w.append(w)

    if non_zero_w:
        mean_w = sum(non_zero_w) / len(non_zero_w)
    else:
        mean_w = 0.0

    # second pass: codons with W_i == 0 get mean(w)
    for codon, W_i in W_by_codon.items():
        if W_i == 0:
            w_by_codon[codon] = mean_w

    return W_by_codon, w_by_codon

