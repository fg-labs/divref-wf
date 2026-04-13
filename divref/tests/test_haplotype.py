"""Tests for shared Hail utilities in haplotype.py."""

import pytest

from divref.haplotype import get_haplo_sequence


def test_get_haplo_sequence_empty_list_raises() -> None:
    """get_haplo_sequence should raise ValueError when given an empty list."""
    with pytest.raises(ValueError, match="at least one variant"):
        get_haplo_sequence(context_size=2, variants=[])


def test_get_haplo_sequence_empty_tuple_raises() -> None:
    """get_haplo_sequence should raise ValueError when given an empty tuple."""
    with pytest.raises(ValueError, match="at least one variant"):
        get_haplo_sequence(context_size=2, variants=())
