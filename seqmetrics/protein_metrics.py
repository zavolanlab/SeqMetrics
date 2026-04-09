from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def translate_cds(seq: str) -> str:
    """Translate nucleotide CDS to protein (no introns, standard code)."""
    dna = Seq(seq.upper())
    protein = dna.translate(to_stop=True)
    return str(protein)

def protein_basic_properties(protein_seq: str):
    """
    Returns:
    - mw: molecular weight
    - aa_comp: dict of amino acid composition (fraction)
    - pI: isoelectric point
    - gravy: hydropathicity
    - aromaticity: fraction aromatic residues
    """
    if not protein_seq:
        return float("nan"), {}, float("nan"), float("nan"), float("nan")

    analysis = ProteinAnalysis(protein_seq)
    mw = analysis.molecular_weight()
    aa_comp = analysis.get_amino_acids_percent()
    pI = analysis.isoelectric_point()
    gravy = analysis.gravy()
    aromaticity = analysis.aromaticity()
    return mw, aa_comp, pI, gravy, aromaticity

def secondary_structure_fraction(protein_seq: str):
    """
    Placeholder: fraction helix, turn, sheet.
    Biopython's ProteinAnalysis.secondary_structure_fraction
    is based on propensity scales, not real structure prediction.
    """
    if not protein_seq:
        return float("nan"), float("nan"), float("nan")
    analysis = ProteinAnalysis(protein_seq)
    helix, turn, sheet = analysis.secondary_structure_fraction()
    return helix, sheet, 1.0 - helix - sheet  # treat remaining as coil
