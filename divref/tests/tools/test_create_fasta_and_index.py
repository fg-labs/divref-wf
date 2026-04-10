"""Tests for get_haplo_sequence in haplotype.py."""

from typing import Any
from unittest.mock import patch

import pytest

hl = pytest.importorskip("hail")

from divref.haplotype import get_haplo_sequence  # noqa: E402


def _create_variant(contig: str, position: int, ref: str, alt: str) -> Any:
    """
    Create a test variant struct in the format expected by get_haplo_sequence.

    Args:
        contig: Chromosome name.
        position: 1-based reference position.
        ref: Reference allele.
        alt: Alternate allele.

    Returns:
        Hail struct with locus and alleles fields.
    """
    return hl.Struct(
        locus=hl.Struct(contig=contig, position=position),
        alleles=[ref, alt],
    )


def _create_reference_mock(reference_sequence: str) -> Any:
    """
    Create a mock for hl.get_sequence that returns substrings of a fixed reference.

    The mock accepts the same arguments as hl.get_sequence and returns the
    appropriate substring of the provided reference string.

    Args:
        reference_sequence: The reference string to use for subsequence extraction.

    Returns:
        A callable that mimics hl.get_sequence using the provided reference.
    """

    def mock_get_sequence(
        _contig: str,
        position: int,
        before: int = 0,
        after: int = 0,
        reference_genome: Any = None,  # noqa: ARG001
    ) -> Any:
        return hl.str(reference_sequence)[position - before : position + after + 1]

    return mock_get_sequence


@pytest.mark.skip(reason="Requires a running Hail/Spark JVM context")
def test_get_haplo_sequence_edge_cases() -> None:
    """Test get_haplo_sequence with SNPs, insertions, and deletions."""
    reference = "01234567891"

    two_snps = [
        _create_variant("chr1", 4, "A", "T"),
        _create_variant("chr1", 6, "G", "C"),
    ]
    two_insertions = [
        _create_variant("chr1", 4, "A", "AT"),
        _create_variant("chr1", 6, "G", "GC"),
    ]
    two_deletions = [
        _create_variant("chr1", 4, "AT", "A"),
        _create_variant("chr1", 7, "GC", "G"),
    ]

    mock_get_sequence = _create_reference_mock(reference)

    with patch("hail.get_sequence", side_effect=mock_get_sequence):
        result_snps = hl.eval(get_haplo_sequence(context_size=2, variants=two_snps))
        assert result_snps == "23T5C78"

        result_insertions = hl.eval(get_haplo_sequence(context_size=2, variants=two_insertions))
        assert result_insertions == "23AT5GC78"

        result_deletions = hl.eval(get_haplo_sequence(context_size=2, variants=two_deletions))
        assert result_deletions == "23A6G91"
