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
    in_gnomad_hgdp_sample_metadata: HailPath = defaults.GNOMAD_HGDP_1KG_SAMPLE_METADATA_HAIL_TABLE,  # noqa: E501
    out_variant_annotation_table: HailPath,
    out_sample_metadata: HailPath,
    locus: str = "chr1:100001-200000",
    sample_proportion: float = 0.01,
    seed: int = 42,
    gcs_credentials_path: Path = Path("~/.config/gcloud/application_default_credentials.json"),
) -> None:
    """
    Extract subsets of gnomAD HGDP/1KG variant annotations and sample metadata for testing.

    Args:
        in_gnomad_hgdp_variant_annotation_table: Path to the gnomAD HGDP/1KG variant annotation
            Hail table.
        in_gnomad_hgdp_sample_metadata: Path to the gnomAD HGDP/1KG sample metadata Hail table.
        out_variant_annotation_table: Output path for the subset variant annotation Hail table.
        out_sample_metadata: Output path for the sample metadata Hail table, stripped to key and
            gnomad_population_inference only.
        locus: Sampling locus for variant filtering.
        sample_proportion: Proportion of variants to sample.
        seed: Random sampling seed.
        gcs_credentials_path: Path to GCS default credentials JSON file.
    """
    if not 0.0 <= sample_proportion <= 1.0:
        raise ValueError("sample_proportion must be between 0.0 and 1.0")

    hail_init(gcs_credentials_path.expanduser())

    va = hl.read_table(in_gnomad_hgdp_variant_annotation_table)
    roi = [hl.parse_locus_interval(locus, reference_genome="GRCh38")]

    logger.info(f"Filtering to {locus} and sampling {sample_proportion} of variants.")
    va_subset = hl.filter_intervals(va, roi).sample(sample_proportion, seed)

    logger.info(f"Writing {va_subset.count()} variants to {out_variant_annotation_table}.")
    va_subset.write(out_variant_annotation_table, overwrite=True)

    sa = hl.read_table(in_gnomad_hgdp_sample_metadata)
    sa = sa.select("gnomad_population_inference").select_globals()

    logger.info(f"Writing {sa.count()} samples to {out_sample_metadata}.")
    sa.naive_coalesce(1).write(out_sample_metadata, overwrite=True)
