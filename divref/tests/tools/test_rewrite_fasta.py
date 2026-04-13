"""Tests for the rewrite_fasta tool."""

from pathlib import Path

from divref.tools.rewrite_fasta import rewrite_fasta


def _write_fasta(path: Path, contigs: list[tuple[str, str]]) -> None:
    with open(path, "w") as f:
        for name, seq in contigs:
            f.write(f">{name}\n{seq}\n")


def test_keeps_canonical_autosomes(tmp_path: Path) -> None:
    _write_fasta(tmp_path / "in.fa", [("chr1", "ACGT"), ("chr22", "TTTT")])
    rewrite_fasta(fasta_path=tmp_path / "in.fa", output_path=tmp_path / "out.fa")
    assert (tmp_path / "out.fa").read_text() == ">chr1\nACGT\n>chr22\nTTTT\n"


def test_keeps_sex_and_mt_chromosomes(tmp_path: Path) -> None:
    _write_fasta(tmp_path / "in.fa", [("chrX", "AAAA"), ("chrY", "CCCC"), ("chrMT", "GGGG")])
    rewrite_fasta(fasta_path=tmp_path / "in.fa", output_path=tmp_path / "out.fa")
    assert (tmp_path / "out.fa").read_text() == ">chrX\nAAAA\n>chrY\nCCCC\n>chrMT\nGGGG\n"


def test_filters_non_canonical(tmp_path: Path) -> None:
    _write_fasta(
        tmp_path / "in.fa",
        [("chr1", "ACGT"), ("chr1_alt", "AAAA"), ("chrUn_gl000220", "TTTT")],
    )
    rewrite_fasta(fasta_path=tmp_path / "in.fa", output_path=tmp_path / "out.fa")
    assert (tmp_path / "out.fa").read_text() == ">chr1\nACGT\n"


def test_mixed_canonical_and_non_canonical(tmp_path: Path) -> None:
    _write_fasta(
        tmp_path / "in.fa",
        [("chr1", "ACGT"), ("chr1_alt", "AAAA"), ("chr2", "TTTT"), ("chrEBV", "CCCC")],
    )
    rewrite_fasta(fasta_path=tmp_path / "in.fa", output_path=tmp_path / "out.fa")
    assert (tmp_path / "out.fa").read_text() == ">chr1\nACGT\n>chr2\nTTTT\n"


def test_empty_fasta(tmp_path: Path) -> None:
    (tmp_path / "in.fa").write_text("")
    rewrite_fasta(fasta_path=tmp_path / "in.fa", output_path=tmp_path / "out.fa")
    assert (tmp_path / "out.fa").read_text() == ""


def test_multiline_sequence_preserved(tmp_path: Path) -> None:
    (tmp_path / "in.fa").write_text(">chr1\nACGT\nTGCA\n>chr1_alt\nAAAA\n")
    rewrite_fasta(fasta_path=tmp_path / "in.fa", output_path=tmp_path / "out.fa")
    assert (tmp_path / "out.fa").read_text() == ">chr1\nACGT\nTGCA\n"
