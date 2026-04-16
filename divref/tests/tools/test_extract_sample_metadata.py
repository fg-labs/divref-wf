"""Tests for the extract_sample_metadata tool."""

from pathlib import Path
from unittest.mock import patch

import hail as hl

from divref.tools.extract_sample_metadata import extract_sample_metadata


def test_extract_sample_metadata(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
) -> None:
    """Happy-path: extract AFs from local test variant and sample data."""
    in_samples = str(datadir / "hgdp_1kg_sample_metadata.ht")
    out_sa = str(tmp_path / "sa.ht")

    with patch("divref.tools.extract_sample_metadata.hail_init"):
        extract_sample_metadata(
            in_gnomad_hgdp_sample_data=in_samples,
            out_sample_metadata=out_sa,
        )

    sa = hl.read_table(out_sa)
    sa_count = sa.count()
    assert sa_count == 4151

    # Each sample should have a pop field
    first_sample = sa.head(1).collect()[0]
    assert hasattr(first_sample, "pop")
    assert isinstance(first_sample.pop, str)
