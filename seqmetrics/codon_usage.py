# codon_usage.py
from __future__ import annotations
from dataclasses import dataclass
from collections import Counter
from typing import Dict, List, Optional, Set
import math
import csv


@dataclass
class CodonUsageTable:
    """
    Stores codon usage counts (e.g. per tissue) and provides
    RSCU, optimal codons, CAI, fraction of optimal codons.
    """
    name: str                    # e.g. "brain" or "liver"
    codon_counts: Counter        # codon -> count
    syn_codons_by_aa: Dict[str, List[str]]

    @classmethod
    def from_row(
        cls,
        name: str,
        codon_to_count: Dict[str, int],
        syn_codons_by_aa: Dict[str, List[str]],
    ) -> "CodonUsageTable":
        return cls(name=name,
                   codon_counts=Counter({c.upper(): int(v) for c, v in codon_to_count.items()}),
                   syn_codons_by_aa=syn_codons_by_aa)

    @staticmethod
    def load_from_table(
        path: str,
        syn_codons_by_aa: Dict[str, List[str]],
        tissue_col: str = "Tissue",
    ) -> Dict[str, "CodonUsageTable"]:
        """
        Load tissue-specific codon usage from a tab-delimited table.

        Expected format:

            Tissue\tAAA\tAAC\tAAG\t...\tTTT\t[other columns...]
            brain \t10 \t5  \t7  \t...\t3
            liver \t20 \t7  \t9  \t...\t1

        - One row per tissue.
        - Column named `tissue_col` (default 'Tissue') holds tissue names.
        - Any header column that is a 3-letter string of A/C/G/T is treated as a codon.
        - Other columns are ignored.
        """
        tables: Dict[str, CodonUsageTable] = {}

        with open(path, "r", encoding="utf-8") as fh:
            reader = csv.reader(fh, delimiter="\t")

            try:
                header_cols = next(reader)
            except StopIteration:
                raise ValueError("Codon-usage file is empty")

            header_cols = [h.strip() for h in header_cols]

            # Tissue column index
            if tissue_col in header_cols:
                tissue_idx = header_cols.index(tissue_col)
            else:
                # fallback: first column
                tissue_idx = 0

            # Codon columns: 3-letter A/C/G/T headers
            codon_indices = []
            codon_names = []
            for i, name in enumerate(header_cols):
                name_up = name.strip().upper()
                if len(name_up) == 3 and all(b in "ACGT" for b in name_up):
                    codon_indices.append(i)
                    codon_names.append(name_up)

            if not codon_indices:
                raise ValueError(
                    "No codon columns found (headers like 'AAA', 'AAC', ...) in "
                    f"{path}"
                )

            # Data rows
            for row in reader:
                if not row or all(not x.strip() for x in row):
                    continue

                # pad short rows
                if len(row) < len(header_cols):
                    row = row + [""] * (len(header_cols) - len(row))

                tissue_name = row[tissue_idx].strip()
                if not tissue_name:
                    continue

                codon_to_count: Dict[str, int] = {}
                for idx, codon in zip(codon_indices, codon_names):
                    if idx >= len(row):
                        continue
                    val = row[idx].strip()
                    if not val:
                        count = 0
                    else:
                        val_clean = val.rstrip(",")
                        try:
                            count = int(val_clean)
                        except ValueError:
                            count = 0
                    codon_to_count[codon] = count

                tables[tissue_name] = CodonUsageTable.from_row(
                    tissue_name,
                    codon_to_count,
                    syn_codons_by_aa,
                )

        return tables

    def rscu(self) -> Dict[str, float]:
        """Relative synonymous codon usage per codon."""
        r: Dict[str, float] = {}
        for aa, codons in self.syn_codons_by_aa.items():
            total = sum(self.codon_counts.get(c, 0) for c in codons)
            k = len(codons)
            if total == 0:
                for c in codons:
                    r[c] = 0.0
            else:
                for c in codons:
                    r[c] = (self.codon_counts.get(c, 0) * k) / total
        return r

    def optimal_codons(self, rscu_threshold: float = 0.9) -> Set[str]:
        """
        For each amino acid, codons with RSCU >= threshold * max_RSCU
        among its synonymous set are considered optimal.
        """
        r = self.rscu()
        opt: Set[str] = set()
        for aa, codons in self.syn_codons_by_aa.items():
            vals = [r.get(c, 0.0) for c in codons]
            if not vals:
                continue
            max_r = max(vals)
            if max_r == 0:
                continue
            for c in codons:
                if r.get(c, 0.0) >= rscu_threshold * max_r:
                    opt.add(c)
        return opt

    def cai(
        self,
        codons_in_seq: List[str],
        min_weight: float = 0.01,
    ) -> float:
        """
        Codon Adaptation Index for a CDS (given as codon list) against this table.
        """
        r = self.rscu()
        w: Dict[str, float] = {}
        for aa, codons in self.syn_codons_by_aa.items():
            vals = [r.get(c, 0.0) for c in codons]
            if not vals:
                continue
            max_r = max(vals)
            if max_r == 0:
                for c in codons:
                    w[c] = min_weight
                continue
            for c in codons:
                w[c] = r.get(c, 0.0) / max_r

        if not codons_in_seq:
            return float("nan")

        logs = []
        for c in codons_in_seq:
            wi = w.get(c.upper(), min_weight)
            if wi <= 0:
                wi = min_weight
            logs.append(math.log(wi))
        return math.exp(sum(logs) / len(logs))

    def fraction_optimal(
        self,
        codons_in_seq: List[str],
        optimal: Optional[Set[str]] = None,
    ) -> float:
        """Fraction of codons in the CDS that are optimal for this table."""
        if not codons_in_seq:
            return 0.0
        if optimal is None:
            optimal = self.optimal_codons()
        opt_set = {c.upper() for c in optimal}
        opt_count = sum(1 for c in codons_in_seq if c.upper() in opt_set)
        return opt_count / len(codons_in_seq)
