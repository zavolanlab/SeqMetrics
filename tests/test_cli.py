from pathlib import Path

from seqmetrics.cli import main


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    lines = []
    for seq_id, seq in records:
        lines.append(f">{seq_id}")
        lines.append(seq)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_cli_writes_tsv(tmp_path, capsys):
    fasta_path = tmp_path / "input.fa"
    write_fasta(fasta_path, [("seq1", "ATGGCTTAA")])

    exit_code = main([str(fasta_path)])

    captured = capsys.readouterr()
    lines = captured.out.strip().splitlines()

    assert exit_code == 0
    assert lines[0].startswith("seq_id\tdescription\tlength_nt")
    assert "seq1" in lines[1]
    assert "aa_M" in lines[0]


def test_cli_creates_output_directory(tmp_path):
    fasta_path = tmp_path / "input.fa"
    output_path = tmp_path / "nested" / "results.tsv"
    write_fasta(fasta_path, [("seq1", "ATGGCTTAA")])

    exit_code = main([str(fasta_path), "--output", str(output_path)])

    assert exit_code == 0
    assert output_path.exists()