"""Tests for the compute_haplotypes tool."""

from pathlib import Path
from unittest.mock import patch

import hail as hl

from divref import defaults
from divref.tools.compute_haplotypes import compute_haplotypes
from divref.tools.extract_gnomad_afs import extract_gnomad_afs
from divref.tools.extract_sample_metadata import extract_sample_metadata


def test_compute_haplotypes(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
) -> None:
    """Happy-path: compute haplotypes from test VCF with gnomAD annotations."""
    # --- setup: produce VA and SA tables from test data ---
    in_sites = str(datadir / "chr1_100001_200000.ht")
    in_samples = str(datadir / "hgdp_1kg_sample_metadata.ht")
    va_path = str(tmp_path / "va.ht")
    sa_path = str(tmp_path / "sa.ht")

    with patch("divref.tools.extract_gnomad_afs.hail_init"):
        extract_gnomad_afs(
            in_gnomad_sites_table=in_sites,
            out_variant_annotation_table=va_path,
            contig="chr1",
            populations=defaults.POPULATIONS,
            reference_genome=defaults.REFERENCE_GENOME,
        )

    with patch("divref.tools.extract_sample_metadata.hail_init"):
        extract_sample_metadata(
            in_gnomad_hgdp_sample_data=in_samples,
            out_sample_metadata=sa_path,
        )

    # --- act ---
    vcf_path = str(datadir / "chr1_100001_200000.vcf.gz")
    output_base = str(tmp_path / "haplos")

    with patch("divref.tools.compute_haplotypes.hl.init"):
        compute_haplotypes(
            vcfs_path=vcf_path,
            gnomad_va_file=va_path,
            gnomad_sa_file=sa_path,
            window_size=5000,
            freq_threshold=0.005,
            output_base=output_base,
            temp_dir=str(tmp_path / "hail_tmp"),
        )

    # --- assert ---
    result = hl.read_table(f"{output_base}.ht")
    result_count = result.count()
    assert result_count == 295

    results: list[hl.Struct] = result.collect()

    # Each haplotype has >= 2 variants (single-variant haplotypes are filtered)
    assert all([len(r.haplotype) >= 2 for r in results])
    assert all([len(r.variants) == len(r.haplotype) for r in results])
    assert all([len(r.gnomad_freqs) == len(r.haplotype) for r in results])

    # Population frequency summary fields are present
    assert all([hasattr(r, "max_pop") for r in results])
    assert all([hasattr(r, "max_empirical_AF") for r in results])
    assert all([hasattr(r, "max_empirical_AC") for r in results])
    assert all([hasattr(r, "all_pop_freqs") for r in results])
    assert all([len(r.all_pop_freqs) > 0 for r in results])
