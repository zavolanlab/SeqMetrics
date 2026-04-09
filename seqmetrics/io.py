# seqmetrics/io.py
from typing import List
from Bio import SeqIO

def read_fasta_cds(path: str):
    """Yield (id, description, cds_seq) for each record in FASTA."""
    for rec in SeqIO.parse(path, "fasta"):
        yield rec.id, rec.description, str(rec.seq).upper()
