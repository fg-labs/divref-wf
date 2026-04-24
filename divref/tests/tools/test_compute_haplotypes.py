"""Tests for the compute_haplotypes tool."""

from pathlib import Path
from unittest.mock import patch

import hail as hl
import pytest

from divref.tools.compute_haplotypes import compute_haplotypes


@pytest.mark.parametrize(
    "variant_freq_threshold,haplotype_freq_threshold,expected_count",
    [
        (0.0, 0.0, 517),
        (0.005, 0.0, 295),
        (0.0, 0.005, 30),
        (0.005, 0.005, 33),
    ],
)
def test_compute_haplotypes(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
    variant_freq_threshold: float,
    haplotype_freq_threshold: float,
    expected_count: int,
) -> None:
    """Happy-path: compute haplotypes from test VCF with gnomAD annotations."""
    # --- act ---
    in_sites = datadir / "chr1_100001_200000.gnomad_afs.ht"
    in_samples = datadir / "hgdp_1kg_sample_metadata.extract.ht"
    vcf_path = datadir / "chr1_100001_200000.vcf.gz"
    output_base = tmp_path / "haplos"

    with patch("divref.tools.compute_haplotypes.hl.init"):
        compute_haplotypes(
            vcfs_path=vcf_path,
            gnomad_va_file=in_sites,
            gnomad_sa_file=in_samples,
            window_size=5000,
            variant_freq_threshold=variant_freq_threshold,
            haplotype_freq_threshold=haplotype_freq_threshold,
            output_base=output_base,
            temp_dir=tmp_path / "hail_tmp",
        )

    # --- assert ---
    result = hl.read_table(f"{output_base}.ht")
    result_count = result.count()
    assert result_count == expected_count

    results: list[hl.Struct] = result.collect()

    # Each haplotype has >= 2 variants (single-variant haplotypes are filtered)
    assert all(len(r.haplotype) >= 2 for r in results)
    assert all(len(r.variants) == len(r.haplotype) for r in results)
    assert all(len(r.gnomad_freqs) == len(r.haplotype) for r in results)

    # Population frequency summary fields are present
    assert all(hasattr(r, "max_pop") for r in results)
    assert all(hasattr(r, "max_empirical_AF") for r in results)
    assert all(hasattr(r, "max_empirical_AC") for r in results)
    assert all(hasattr(r, "all_pop_freqs") for r in results)
    assert all(len(r.all_pop_freqs) > 0 for r in results)


def test_compute_haplotypes_no_variants(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
) -> None:
    """All variants are filtered out."""
    # --- act ---
    in_sites = datadir / "chr1_100001_200000.gnomad_afs.ht"
    in_samples = datadir / "hgdp_1kg_sample_metadata.extract.ht"
    vcf_path = datadir / "chr1_100001_200000.vcf.gz"
    output_base = tmp_path / "haplos"

    with (
        patch("divref.tools.compute_haplotypes.hl.init"),
        pytest.raises(ValueError, match="No variants found with minimum population AF"),
    ):
        compute_haplotypes(
            vcfs_path=vcf_path,
            gnomad_va_file=in_sites,
            gnomad_sa_file=in_samples,
            window_size=5000,
            variant_freq_threshold=1,
            haplotype_freq_threshold=0,
            output_base=output_base,
            temp_dir=tmp_path / "hail_tmp",
        )
