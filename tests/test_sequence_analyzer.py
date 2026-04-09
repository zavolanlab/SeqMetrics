import math
from pathlib import Path

from seqmetrics import CodonUsageTable, SequenceAnalyzer
from seqmetrics.codon_definitions import SYN_CODONS_BY_AA


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    lines = []
    for seq_id, seq in records:
        lines.append(f">{seq_id}")
        lines.append(seq)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_analyze_fasta_without_optional_tables(tmp_path):
    fasta_path = tmp_path / "input.fa"
    write_fasta(fasta_path, [("seq1", "ATGGCTTAA")])

    analyzer = SequenceAnalyzer()
    rows = analyzer.analyze_fasta(str(fasta_path))

    assert len(rows) == 1
    row = rows[0]
    assert row["seq_id"] == "seq1"
    assert row["length_nt"] == 9
    assert row["length_aa"] == 2
    assert row["tissue"] is None
    assert row["species"] is None
    assert math.isnan(row["tissue_cai"])
    assert math.isnan(row["species_tai"])
    assert row["aa_A"] == 50.0
    assert row["aa_M"] == 50.0


def test_analyze_fasta_with_tissue_usage(tmp_path):
    fasta_path = tmp_path / "input.fa"
    table_path = tmp_path / "tissue_usage.tsv"
    write_fasta(fasta_path, [("seq1", "AAAGCT")])
    table_path.write_text(
        "Tissue\tAAA\tAAG\tGCT\tGCC\n"
        "brain\t100\t0\t80\t20\n",
        encoding="utf-8",
    )

    tissue_usage = CodonUsageTable.load_from_table(
        str(table_path),
        syn_codons_by_aa=SYN_CODONS_BY_AA,
    )
    analyzer = SequenceAnalyzer(tissue_usage=tissue_usage, default_tissue="brain")

    rows = analyzer.analyze_fasta(str(fasta_path))

    assert len(rows) == 1
    row = rows[0]
    assert row["tissue"] == "brain"
    assert row["tissue_cai"] == 1.0
    assert row["frac_opt_codons"] == 1.0