"""Tool to extract gnomAD sample metadata for the DivRef pipeline."""

from pathlib import Path

import hail as hl
from fgpyo.io import assert_path_is_writable

from divref.alias import HailPath
from divref.hail import hail_init


def extract_sample_metadata(
    *,
    in_gnomad_hgdp_sample_data: HailPath,
    out_sample_metadata: Path,
    gcs_credentials_path: Path = Path("~/.config/gcloud/application_default_credentials.json"),
) -> None:
    """
    Extract sample metadata for downstream pipeline tools.

    Args:
        in_gnomad_hgdp_sample_data: Path to the gnomAD HGDP/1KG sample metadata table.
        out_sample_metadata: Output path for the sample metadata Hail table.
        gcs_credentials_path: Path to GCS default credentials JSON file.
    """
    assert_path_is_writable(out_sample_metadata)

    hail_init(gcs_credentials_path.expanduser())

    sa = hl.read_table(in_gnomad_hgdp_sample_data).select_globals()
    sa = sa.select(pop=sa.gnomad_population_inference.pop)
    sa.naive_coalesce(4).write(str(out_sample_metadata), overwrite=True)
