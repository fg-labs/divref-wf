"""Tests for the compute_haplotype_statistics tool."""

from pathlib import Path
from unittest.mock import patch

import hail as hl

from divref.tools.compute_haplotype_statistics import compute_haplotype_statistics


def test_compute_haplotype_statistics(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
) -> None:
    """Happy-path: compute haplotype and gnomAD statistics and verify output TSV structure."""
    in_haplotypes = str(datadir / "chr1_100001_200000_haplotypes.ht")
    in_sites = str(datadir / "chr1_100001_200000.gnomad_afs.ht")
    output_base = tmp_path / "stats"

    with patch("divref.tools.compute_haplotype_statistics.hl.init"):
        compute_haplotype_statistics(
            haplotypes_table_path=in_haplotypes,
            gnomad_va_file=in_sites,
            window_sizes=[5000],
            frequency_cutoffs=[0.005],
            output_base=output_base,
        )

    # --- assert HGDP output ---
    hgdp_tsv = output_base.with_suffix(".hgdp.tsv")
    assert hgdp_tsv.exists()

    hgdp_lines = hgdp_tsv.read_text().splitlines()
    assert hgdp_lines[0] == "frequency\twindow_size\thgdp_haplotype_count"

    # One row per (frequency, window_size) combination
    assert len(hgdp_lines) == 2  # header + 1 data row (1 freq × 1 window)

    freq, window, count = hgdp_lines[1].split("\t")
    assert float(freq) == 0.005
    assert int(window) == 5000
    assert int(count) == 24

    # --- assert gnomAD output ---
    gnomad_tsv = output_base.with_suffix(".gnomad.tsv")
    assert gnomad_tsv.exists()

    gnomad_lines = gnomad_tsv.read_text().splitlines()
    assert gnomad_lines[0] == "frequency\tgnomad_variant_count"

    # One row per frequency cutoff
    assert len(gnomad_lines) == 2  # header + 1 data row (1 freq)

    freq, gnomad_count = gnomad_lines[1].split("\t")
    assert float(freq) == 0.005
    assert int(gnomad_count) == 473
