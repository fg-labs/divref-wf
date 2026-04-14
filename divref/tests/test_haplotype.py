"""Tests for shared Hail utilities in haplotype.py."""

from pathlib import Path
from typing import Any

import hail as hl
import pytest

from divref.haplotype import get_haplo_sequence
from divref.haplotype import split_haplotypes
from divref.haplotype import to_hashable_items
from divref.haplotype import variant_distance

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def _make_variant(position: int, ref: str, alt: str) -> Any:
    return hl.Struct(locus=hl.Struct(contig="chr1", position=position), alleles=[ref, alt])


def _make_haplotype_table(variant_positions: list[tuple[str, int, str, str]]) -> Any:
    variant_type = hl.tstruct(
        locus=hl.tstruct(contig=hl.tstr, position=hl.tint32), alleles=hl.tarray(hl.tstr)
    )
    row_type = hl.tstruct(
        variants=hl.tarray(variant_type),
        haplotype=hl.tarray(hl.tstr),
        gnomad_freqs=hl.tarray(hl.tfloat64),
    )
    variants = [
        {"locus": {"contig": contig, "position": pos}, "alleles": [ref, alt]}
        for contig, pos, ref, alt in variant_positions
    ]
    return hl.Table.parallelize(
        [
            {
                "variants": variants,
                "haplotype": [str(i) for i in range(len(variants))],
                "gnomad_freqs": [0.1] * len(variants),
            }
        ],
        schema=row_type,
    )


# ---------------------------------------------------------------------------
# get_haplo_sequence
# ---------------------------------------------------------------------------


def test_get_haplo_sequence_single(
    datadir: Path,
    hail_reference_genome: Any,
    hail_context: None,  # noqa: ARG001
) -> None:
    """get_haplo_sequence should return the correct haplotype sequence."""
    test_fasta: Path = datadir / "test.fa"
    test_fai: Path = datadir / "test.fa.fai"

    hail_reference_genome.add_sequence(str(test_fasta), str(test_fai))

    variant = _make_variant(position=100, ref="A", alt="C")
    haplo_seq = get_haplo_sequence(
        context_size=2, variants=[variant], reference_genome=hail_reference_genome.name
    )
    assert hl.eval(haplo_seq) == "CCCTC"


def test_get_haplo_sequence_empty_list_raises() -> None:
    """get_haplo_sequence should raise ValueError when given an empty list."""
    with pytest.raises(ValueError, match="at least one variant"):
        get_haplo_sequence(context_size=2, variants=[])


def test_get_haplo_sequence_empty_tuple_raises() -> None:
    """get_haplo_sequence should raise ValueError when given an empty tuple."""
    with pytest.raises(ValueError, match="at least one variant"):
        get_haplo_sequence(context_size=2, variants=())


# ---------------------------------------------------------------------------
# to_hashable_items
# ---------------------------------------------------------------------------


def test_to_hashable_items_empty() -> None:
    assert to_hashable_items({}) == ()


def test_to_hashable_items_single_entry() -> None:
    assert to_hashable_items({"key": "value"}) == (("key", "value"),)


def test_to_hashable_items_sorted_by_key() -> None:
    assert to_hashable_items({"b": 2, "a": 1, "c": 3}) == (("a", 1), ("b", 2), ("c", 3))


# ---------------------------------------------------------------------------
# variant_distance
# ---------------------------------------------------------------------------


def test_variant_distance_adjacent_snps(hail_context: None) -> None:  # noqa: ARG001
    # SNP at 100, next SNP at 101: distance = 101 - 100 - len("A") = 0
    assert (
        hl.eval(variant_distance(_make_variant(100, "A", "T"), _make_variant(101, "C", "G"))) == 0
    )


def test_variant_distance_snps_with_gap(hail_context: None) -> None:  # noqa: ARG001
    # SNP at 100, next SNP at 103: 2 reference bases separate them
    assert (
        hl.eval(variant_distance(_make_variant(100, "A", "T"), _make_variant(103, "C", "G"))) == 2
    )


def test_variant_distance_deletion_closes_gap(hail_context: None) -> None:  # noqa: ARG001
    # Deletion AT→A at 100 (consumes 2 ref bases), next variant at 102: distance = 0
    assert (
        hl.eval(variant_distance(_make_variant(100, "AT", "A"), _make_variant(102, "C", "G"))) == 0
    )


# ---------------------------------------------------------------------------
# split_haplotypes
# ---------------------------------------------------------------------------


def test_split_haplotypes_no_split_needed(hail_context: None) -> None:  # noqa: ARG001
    # All variants within window_size=200; haplotype is kept intact as one row
    ht = _make_haplotype_table([
        ("chr1", 100, "A", "T"),
        ("chr1", 150, "C", "G"),
        ("chr1", 190, "G", "A"),
    ])
    rows = split_haplotypes(ht, window_size=200).collect()
    assert len(rows) == 1
    assert len(rows[0].variants) == 3


def test_split_haplotypes_splits_at_large_gap(hail_context: None) -> None:  # noqa: ARG001
    # Gap between positions 101 and 500 (398 bases) exceeds window_size=200;
    # results in two sub-haplotypes: [v0, v1] and [v2, v3]
    ht = _make_haplotype_table([
        ("chr1", 100, "A", "T"),
        ("chr1", 101, "C", "G"),
        ("chr1", 500, "G", "A"),
        ("chr1", 501, "T", "C"),
    ])
    rows = sorted(
        split_haplotypes(ht, window_size=200).collect(),
        key=lambda r: r.variants[0].locus.position,
    )
    assert len(rows) == 2
    assert [v.locus.position for v in rows[0].variants] == [100, 101]
    assert [v.locus.position for v in rows[1].variants] == [500, 501]


def test_split_haplotypes_discards_singleton_segment(hail_context: None) -> None:  # noqa: ARG001
    # Gap after position 100 isolates it as a singleton (discarded);
    # only the two-variant segment [500, 501] is kept
    ht = _make_haplotype_table([
        ("chr1", 100, "A", "T"),
        ("chr1", 500, "C", "G"),
        ("chr1", 501, "G", "A"),
    ])
    rows = split_haplotypes(ht, window_size=200).collect()
    assert len(rows) == 1
    assert [v.locus.position for v in rows[0].variants] == [500, 501]
