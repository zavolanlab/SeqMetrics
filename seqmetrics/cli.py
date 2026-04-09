# seqmetrics/cli.py

import argparse
import sys
import csv
from pathlib import Path

from .sequence_analyzer import SequenceAnalyzer
from .codon_usage import CodonUsageTable
from .codon_definitions import SYN_CODONS_BY_AA
from .trna_weights import TRNAWeightTable

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compute sequence metrics from CDS FASTA."
    )
    parser.add_argument("fasta", help="Input FASTA file of CDS")
    parser.add_argument(
        "-o",
        "--output",
        help="Output TSV file (default: stdout)",
        default="-",
    )
    parser.add_argument(
        "-s",
        "--species",
        help="Species name (used only if a tRNA table is configured)",
        default=None,
    )
    parser.add_argument(
        "-t",
        "--tissue",
        help="Tissue name (must match a tissue in the codon-usage file)",
        default=None,
    )
    parser.add_argument(
        "--tissue-usage-file",
        help="Path to tissue codon-usage table (header: Tissue <codon...>)",
        default=None,
    )
    parser.add_argument(
        "--trna-counts",
        help="File of tRNA gene copy number from gtRNAdb",
        default=None,
    )

    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    # 1) Load tissue codon usage if a file is provided
    if args.tissue_usage_file is not None:
        tissue_usage = CodonUsageTable.load_from_table(
            args.tissue_usage_file,
            syn_codons_by_aa=SYN_CODONS_BY_AA,
        )
    else:
        tissue_usage = None  # no tissue-specific metrics

    # 2) Load tRNA table here if you want species_tai filled.
    if args.trna_counts is not None:
        trna_table = TRNAWeightTable.from_trna_gene_file(
            species=args.species or "H_sapiens",
            path=args.trna_counts,
            anticodon_col="AntiCodon",
            aa_col="Type",      
            note_col="Note",
            high_conf_phrase="high confidence set",
        )
    else:
        trna_table = None

    # 3) Construct analyzer
    sa = SequenceAnalyzer(
        tissue_usage=tissue_usage,
        trna_table=trna_table,
        default_tissue=args.tissue,
        default_species=args.species,
    )

    # 4) Analyze FASTA
    rows = sa.analyze_fasta(args.fasta, tissue=args.tissue)
    if not rows:
        return 0

    fieldnames = list(rows[0].keys())

    if args.output == "-" or args.output is None:
        out_fh = sys.stdout
    else:
        output_path = Path(args.output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        out_fh = output_path.open("w", newline="", encoding="utf-8")

    try:
        writer = csv.DictWriter(out_fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    finally:
        if out_fh is not sys.stdout:
            out_fh.close()

    return 0


if __name__ == "__main__":
    sys.exit(main())
