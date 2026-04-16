"""Tests for the create_gnomad_sites_vcf tool."""

from pathlib import Path
from unittest.mock import patch

import hail as hl

from divref import defaults
from divref.tools.create_gnomad_sites_vcf import create_gnomad_sites_vcf
from divref.tools.extract_gnomad_afs import extract_gnomad_afs


def test_create_gnomad_sites_vcf(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
) -> None:
    """Happy-path: filter extracted AFs and export as VCF."""
    # --- setup: run extract_gnomad_afs to produce the input table ---
    in_sites = str(datadir / "chr1_100001_200000.ht")
    extracted_table = str(tmp_path / "extracted.ht")

    with patch("divref.tools.extract_gnomad_afs.hail_init"):
        extract_gnomad_afs(
            in_gnomad_sites_table=in_sites,
            out_variant_annotation_table=extracted_table,
            contig="chr1",
            populations=defaults.POPULATIONS,
            reference_genome=defaults.REFERENCE_GENOME,
        )

    # --- act: run create_gnomad_sites_vcf ---
    output_vcf = str(tmp_path / "output.vcf.bgz")

    with patch("divref.tools.create_gnomad_sites_vcf.hail_init"):
        create_gnomad_sites_vcf(
            sites_table_path=extracted_table,
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
    extracted = hl.read_table(extracted_table)
    extracted_count = extracted.count()
    assert vcf_count <= extracted_count
