# sequence_analyzer.py

from __future__ import annotations
from dataclasses import asdict, dataclass
from typing import Dict, Any, List, Optional

from .sequences import SequenceSet, CDSRecord
from .codon_usage import CodonUsageTable
from .trna_weights import TRNAWeightTable
from .codon_definitions import SYN_CODONS_BY_AA 


@dataclass
class SequenceMetrics:
    seq_id: str
    description: str
    length_nt: int
    A: float
    C: float
    G: float
    T: float

    tissue: Optional[str]
    tissue_cai: float
    frac_opt_codons: float

    species: Optional[str]
    species_tai: float

    length_aa: int
    mw: float
    pI: float
    gravy: float
    aromaticity: float
    helix_frac: float
    sheet_frac: float
    coil_frac: float

    aa_composition: Dict[str, float]

class SequenceAnalyzer:
    """
    Connects:
      - sequences (CDSRecord / SequenceSet)
      - tissue-specific codon usage (CodonUsageTable)
      - species-specific tRNA weights (TRNAWeightTable)
      - protein properties
    """

    def __init__(
        self,
        tissue_usage: Optional[Dict[str, CodonUsageTable]] = None,
        trna_table: Optional[TRNAWeightTable] = None,
        default_tissue: Optional[str] = None,
        default_species: Optional[str] = None,
    ):
        """
        Args:
            tissue_usage: mapping tissue -> CodonUsageTable
            trna_table: TRNAWeightTable for one species
            default_tissue: tissue to use if none specified per-call
            default_species: species name to annotate tAI column
        """
        self.tissue_usage = tissue_usage or {}
        self.trna_table = trna_table
        self.default_tissue = default_tissue
        self.default_species = default_species or (
            trna_table.species if trna_table else None
        )

        # cache of optimal codons per tissue
        self._optimal_codons_cache: Dict[str, set] = {}

    # ------------- internal helpers -------------

    def _get_tissue_table(self, tissue: Optional[str]) -> Optional[CodonUsageTable]:
        # If no tissue usage was configured at all, skip tissue-specific metrics
        if not self.tissue_usage:
            return None

        # Prefer explicit tissue argument; fall back to default_tissue
        if tissue is None:
            tissue = self.default_tissue
        if tissue is None:
            return None

        table = self.tissue_usage.get(tissue)
        # If this tissue is not available, quietly skip tissue-specific metrics
        if table is None:
            return None
        return table

    def _get_optimal_codons(self, tissue: str) -> set:
        if tissue not in self._optimal_codons_cache:
            table = self._get_tissue_table(tissue)
            self._optimal_codons_cache[tissue] = table.optimal_codons()
        return self._optimal_codons_cache[tissue]

    def _analyze_record(
        self,
        rec: CDSRecord,
        tissue: Optional[str] = None,
    ) -> SequenceMetrics:
        # basic sequence features
        nt_comp = rec.nucleotide_composition()
        codons = rec.codons()
        length_nt = len(rec.cds)

        # CAI / fraction optimal (tissue-specific)
        tissue_table = self._get_tissue_table(tissue)
        if tissue_table is not None:
            cai_val = tissue_table.cai(codons)
            optimal_codons = self._get_optimal_codons(tissue_table.name)
            frac_opt_val = tissue_table.fraction_optimal(codons, optimal_codons)
            tissue_name = tissue_table.name
        else:
            cai_val = float("nan")
            frac_opt_val = float("nan")
            tissue_name = None

        # tAI (species-specific via existing TRNAWeightTable)
        if self.trna_table is not None:
            tai_val = self.trna_table.tAI(codons)
            species_name = self.trna_table.species
        else:
            tai_val = float("nan")
            species_name = None

        # protein properties
        prot_props = rec.protein_properties()

        return SequenceMetrics(
            seq_id=rec.seq_id,
            description=rec.description,
            length_nt=length_nt,
            A=nt_comp["A"],
            C=nt_comp["C"],
            G=nt_comp["G"],
            T=nt_comp["T"],
            tissue=tissue_name,
            tissue_cai=cai_val,
            frac_opt_codons=frac_opt_val,
            species=species_name,
            species_tai=tai_val,
            length_aa=prot_props["length_aa"],
            mw=prot_props["mw"],
            pI=prot_props["pI"],
            gravy=prot_props["gravy"],
            aromaticity=prot_props["aromaticity"],
            helix_frac=prot_props["helix_frac"],
            sheet_frac=prot_props["sheet_frac"],
            coil_frac=prot_props["coil_frac"],
            aa_composition=prot_props["aa_composition"],
        )

    # ------------- public API -------------

    def analyze_fasta(
        self,
        fasta_path: str,
        tissue: Optional[str] = None,
    ) -> List[Dict[str, Any]]:
        """
        Analyze all CDS in a FASTA file and return a list of dicts
        (one per sequence), ready for conversion to a pandas DataFrame.
        """
        seqs = SequenceSet.from_fasta(fasta_path)
        rows: List[Dict[str, Any]] = []

        for rec in seqs:
            metrics = self._analyze_record(rec, tissue=tissue)
            base = asdict(metrics)

            # flatten aa_composition into columns aa_X
            aa_comp = base.pop("aa_composition", {})
            for aa, frac in aa_comp.items():
                base[f"aa_{aa}"] = frac

            rows.append(base)

        return rows
