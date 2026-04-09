# seqmetrics/sequences.py
from dataclasses import dataclass
from typing import Dict, Iterable, List, Union
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis


@dataclass
class CDSRecord:
    seq_id: str
    description: str
    cds: str

    def codons(self) -> List[str]:
        s = self.cds.upper()
        return [s[i:i+3] for i in range(0, len(s), 3) if len(s[i:i+3]) == 3]

    def nucleotide_composition(self) -> Dict[str, float]:
        s = self.cds.upper().replace("N", "")
        length = len(s) or 1
        return {
            "A": s.count("A") / length,
            "C": s.count("C") / length,
            "G": s.count("G") / length,
            "T": s.count("T") / length,
        }

    def translate(self) -> str:
        dna = Seq(self.cds.upper())
        return str(dna.translate(to_stop=True))

    def protein_properties(self) -> Dict[str, Union[float, Dict[str, float]]]:
        prot = self.translate()
        if not prot:
            return {
                "length_aa": 0,
                "mw": float("nan"),
                "pI": float("nan"),
                "gravy": float("nan"),
                "aromaticity": float("nan"),
                "helix_frac": float("nan"),
                "sheet_frac": float("nan"),
                "coil_frac": float("nan"),
                "aa_composition": {},
            }
        pa = ProteinAnalysis(prot)
        helix, turn, sheet = pa.secondary_structure_fraction()
        return {
            "length_aa": len(prot),
            "mw": pa.molecular_weight(),
            "pI": pa.isoelectric_point(),
            "gravy": pa.gravy(),
            "aromaticity": pa.aromaticity(),
            "helix_frac": helix,
            "sheet_frac": sheet,
            "coil_frac": 1.0 - helix - sheet,
            "aa_composition": pa.amino_acids_percent,  # <- updated line
        }

class SequenceSet:
    def __init__(self, records: Iterable[CDSRecord]):
        self.records: List[CDSRecord] = list(records)

    @classmethod
    def from_fasta(cls, path: str) -> "SequenceSet":
        recs = []
        for rec in SeqIO.parse(path, "fasta"):
            recs.append(CDSRecord(rec.id, rec.description, str(rec.seq).upper()))
        return cls(recs)

    def __iter__(self):
        return iter(self.records)

    def __len__(self):
        return len(self.records)
