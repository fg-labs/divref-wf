"""Tests for the compute_variation_ratios tool."""

from pathlib import Path
from unittest.mock import patch

import hail as hl
import pytest

from divref.tools.compute_variation_ratios import compute_variation_ratios


def test_compute_variation_ratios(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
) -> None:
    """Happy-path: compute per-sample variant site counts and verify output table structure."""
    output_ht = str(tmp_path / "variation_ratios.ht")

    with patch("divref.tools.compute_variation_ratios.hl.init"):
        compute_variation_ratios(
            vcfs_path=str(datadir / "chr1_100001_200000.vcf.gz"),
            gnomad_va_file=str(datadir / "chr1_100001_200000.gnomad_afs.ht"),
            gnomad_sa_file=str(datadir / "hgdp_1kg_sample_metadata.extract.ht"),
            output_ht=output_ht,
            frequency_thresholds=[0.005],
        )

    result = hl.read_table(output_ht)

    # One row per sample in the VCF
    assert result.count() == 4091

    rows: list[hl.Struct] = result.collect()

    # All samples have a population label
    known_pops = {"afr", "amr", "eas", "fin", "mid", "nfe", "oth", "sas"}
    assert all(r.pop in known_pops for r in rows)

    # The counts struct contains the expected per-threshold field
    assert all(hasattr(r.counts, "n_sites_above_0_005") for r in rows)

    # At least some samples have sites above the threshold
    assert any(r.counts.n_sites_above_0_005 > 0 for r in rows)


@pytest.mark.parametrize(
    "thresholds,error",
    [
        ([-0.001], "All frequency_thresholds must be in"),
        ([1.001], "All frequency_thresholds must be in"),
        ([0.005, 0.0050], "produce duplicate output field names"),
    ],
)
def test_compute_variation_ratios_invalid_thresholds(
    thresholds: list[float],
    error: str,
) -> None:
    """Error-path: invalid frequency_thresholds raise ValueError before any Hail work."""
    with pytest.raises(ValueError, match=error):
        compute_variation_ratios(
            vcfs_path="dummy.vcf.gz",
            gnomad_va_file="dummy.ht",
            gnomad_sa_file="dummy.ht",
            output_ht="dummy_out.ht",
            frequency_thresholds=thresholds,
        )
