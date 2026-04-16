"""Tool to extract gnomAD sample metadata for the DivRef pipeline."""

from pathlib import Path

import hail as hl

from divref.alias import HailPath
from divref.hail import hail_init

DEFAULT_POPULATIONS: list[str] = ["afr", "amr", "eas", "sas", "nfe"]
"""Default populations to extract AFs for."""


def extract_sample_metadata(
    *,
    in_gnomad_hgdp_sample_data: HailPath,
    out_sample_metadata: HailPath,
    gcs_credentials_path: Path = Path("~/.config/gcloud/application_default_credentials.json"),
) -> None:
    """
    Extract sample metadata for downstream pipeline tools.

    Args:
        in_gnomad_hgdp_sample_data: Path to the gnomAD HGDP/1KG sample metadata table.
        out_sample_metadata: Output path for the sample metadata Hail table.
        gcs_credentials_path: Path to GCS default credentials JSON file.
    """
    hail_init(gcs_credentials_path.expanduser())

    sa = hl.read_table(in_gnomad_hgdp_sample_data).select_globals()
    sa = sa.select(pop=sa.gnomad_population_inference.pop)
    sa.naive_coalesce(4).write(out_sample_metadata, overwrite=True)
