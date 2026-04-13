"""Shared utilities for Hail-based DivRef pipeline tools."""

from typing import Any
from typing import Hashable
from typing import TypeVar

import hail as hl

HailPath = str
"""Type alias for filesystem paths accepted by Hail: local, GCS (gs://), or HDFS (hdfs://)."""

_V = TypeVar("_V", bound=Hashable)
"""Type variable for hashable dictionary values used in to_hashable_items."""


def to_hashable_items(d: dict[str, _V]) -> tuple[tuple[str, _V], ...]:
    """
    Convert a dictionary to a sorted tuple of items for use as a hashable key.

    Args:
        d: Dictionary with hashable values to convert.

    Returns:
        Sorted tuple of (key, value) pairs.
    """
    return tuple(sorted(d.items()))


def get_haplo_sequence(context_size: int, variants: Any) -> Any:
    """
    Construct a haplotype sequence string with flanking genomic context.

    Builds a sequence by combining alternate alleles from each variant with
    intervening reference sequence, bounded by context_size flanking bases on
    each side.

    Args:
        context_size: Number of reference bases to include flanking each end.
        variants: Hail array expression of variant structs with locus and alleles fields.
            Must contain at least one variant.

    Returns:
        Hail string expression representing the full haplotype sequence.

    Raises:
        ValueError: If variants is a Python sequence with no elements.
    """
    if isinstance(variants, (list, tuple)) and len(variants) == 0:
        raise ValueError(
            "get_haplo_sequence requires at least one variant; received an empty sequence"
        )
    sorted_variants = hl.sorted(variants, key=lambda x: x.locus.position)
    min_variant = sorted_variants[0]
    max_variant = sorted_variants[-1]
    min_pos = min_variant.locus.position
    max_pos = max_variant.locus.position
    max_variant_size = hl.len(max_variant.alleles[0])
    full_context = hl.get_sequence(
        min_variant.locus.contig,
        min_pos,
        before=context_size,
        after=(max_pos - min_pos + max_variant_size + context_size - 1),
        reference_genome="GRCh38",
    )

    # (min_pos - index_translation) equals context_size, mapping locus positions to string indices
    index_translation = min_pos - context_size

    def get_chunk_until_next_variant(i: Any) -> Any:
        v = sorted_variants[i]
        variant_size = hl.len(v.alleles[0])
        reference_buffer_size = hl.if_else(
            i == hl.len(sorted_variants) - 1,
            context_size,
            sorted_variants[i + 1].locus.position - (v.locus.position + variant_size),
        )
        start = v.locus.position - index_translation + variant_size
        return v.alleles[1] + full_context[start : start + reference_buffer_size]

    return full_context[:context_size] + hl.delimit(
        hl.range(hl.len(sorted_variants)).map(get_chunk_until_next_variant),
        "",
    )


def variant_distance(v1: Any, v2: Any) -> Any:
    """
    Calculate the number of reference bases between two variants.

    For example: 1:1:A:T and 1:3:A:T have distance 1 (one base between them).
    1:1:AA:T and 1:3:A:T have distance 0 (deletion closes the gap).

    Args:
        v1: First variant Hail struct with locus and alleles fields.
        v2: Second variant Hail struct with locus and alleles fields.

    Returns:
        Hail int32 expression for the number of reference bases between v1 and v2.
    """
    return v2.locus.position - v1.locus.position - hl.len(v1.alleles[0])


def split_haplotypes(ht: Any, window_size: int) -> Any:
    """
    Split multi-variant haplotypes at gaps of at least `window_size` bases.

    Haplotypes spanning variants further than or equal to `window_size` bases apart are broken
    into sub-haplotypes at those gaps. Sub-haplotypes with fewer than two variants
    are discarded.

    Args:
        ht: Hail table with variants, haplotype, and gnomad_freqs array fields.
        window_size: Maximum reference bases allowed between adjacent variants in
            a haplotype segment.

    Returns:
        Hail table with haplotypes exploded into sub-haplotypes by window.
    """
    breakpoints = hl.range(1, hl.len(ht.variants)).filter(
        lambda i: (i == 0) | (variant_distance(ht.variants[i - 1], ht.variants[i]) >= window_size)
    )

    def get_range(i: Any) -> Any:
        start_index = hl.if_else(i == 0, 0, breakpoints[i - 1])
        end_index = hl.if_else(i == hl.len(breakpoints), hl.len(ht.variants), breakpoints[i])
        return hl.range(start_index, end_index)

    split_hap_indices = (
        hl.range(0, hl.len(breakpoints) + 1).map(get_range).filter(lambda r: hl.len(r) > 1)
    )
    ht = ht.annotate(haplotype_indices=split_hap_indices)
    ht = ht.explode("haplotype_indices")
    ht = ht.annotate(
        haplotype=ht.haplotype_indices.map(lambda i: ht.haplotype[i]),
        variants=ht.haplotype_indices.map(lambda i: ht.variants[i]),
        gnomad_freqs=ht.haplotype_indices.map(lambda i: ht.gnomad_freqs[i]),
    )
    return ht.drop("haplotype_indices")
