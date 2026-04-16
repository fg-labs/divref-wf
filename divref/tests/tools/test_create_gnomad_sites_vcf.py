"""Tests for the create_gnomad_sites_vcf tool."""

from pathlib import Path
from unittest.mock import patch

import hail as hl

from divref import defaults
from divref.tools.create_gnomad_sites_vcf import create_gnomad_sites_vcf


def test_create_gnomad_sites_vcf(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
) -> None:
    """Happy-path: filter extracted AFs and export as VCF."""
    # --- act: run create_gnomad_sites_vcf ---
    in_sites = str(datadir / "chr1_100001_200000.gnomad_afs.ht")
    output_vcf = str(tmp_path / "output.vcf.bgz")

    with patch("divref.tools.create_gnomad_sites_vcf.hail_init"):
        create_gnomad_sites_vcf(
            sites_table_path=in_sites,
            output_vcf_path=output_vcf,
            min_popmax=0.01,
        )

    # --- assert: VCF file is valid ---
    vcf = hl.import_vcf(output_vcf, reference_genome=defaults.REFERENCE_GENOME)
    vcf_count = vcf.count_rows()
    assert vcf_count == 333

    # INFO fields should contain population-prefixed frequency entries
    info_fields = list(vcf.row.info.dtype)
    for pop in defaults.POPULATIONS:
        assert any(f.startswith(f"{pop}_") for f in info_fields), (
            f"Expected INFO fields for population {pop}"
        )

    # VCF should have equal or fewer variants than the unfiltered extracted table
    va = hl.read_table(in_sites)
    assert vcf_count <= va.count()
