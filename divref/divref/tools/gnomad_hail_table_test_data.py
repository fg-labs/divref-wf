"""Tool to fetch a small subset of a gnomAD sites hail table and HGDP sample information."""

import logging
from pathlib import Path

import hail as hl

from divref import defaults
from divref.alias import HailPath
from divref.hail import hail_init

logger = logging.getLogger(__name__)


def gnomad_hail_table_test_data(
    *,
    in_gnomad_hgdp_variant_annotation_table: HailPath = defaults.GNOMAD_HGDP_1KG_VARIANT_ANNOTATION_HAIL_TABLE,  # noqa: E501
    out_variant_annotation_table: HailPath,
    locus: str = "chr1:100001-200000",
    sample_proportion: float = 0.01,
    seed: int = 42,
    gcs_credentials_path: Path = Path("~/.config/gcloud/application_default_credentials.json"),
) -> None:
    """
    Extracts a subset of the gnomAD HGDP/1KG variant annotations for testing.

    Args:
        in_gnomad_hgdp_variant_annotation_table: Path to the gnomAD HGDP/1KG variant annotation
            Hail table.
        in_gnomad_hgdp_sample_metadata: Path to the gnomAD HGDP/1KG sample metadata Hail table.
        out_variant_annotation_table: Output path for the subset variant annotation Hail table.
        out_sample_metadata: Output path for the subset sample metadata Hail table.
        locus: Sampling locus.
        sample_proportion: Proportion of variants to sample.
        seed: Random sampling seed.
        gcs_credentials_path: Path to GCS default credentials JSON file.
    """
    hail_init(gcs_credentials_path.expanduser())

    va = hl.read_table(in_gnomad_hgdp_variant_annotation_table)
    roi = [hl.parse_locus_interval(locus, reference_genome="GRCh38")]

    logger.info(f"Filtering to {locus} and sampling {sample_proportion} of variants.")
    va_subset = hl.filter_intervals(va, roi).sample(sample_proportion, seed)
    logger.info(f"{va_subset.count()} variants found.")

    va_subset.write(out_variant_annotation_table, overwrite=True)
