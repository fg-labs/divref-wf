"""Tests for the extract_gnomad_afs tool."""

from pathlib import Path
from unittest.mock import patch

import hail as hl

from divref.tools.extract_gnomad_afs import DEFAULT_POPULATIONS
from divref.tools.extract_gnomad_afs import extract_gnomad_afs


def test_extract_gnomad_afs(
    hail_context: None,  # noqa: ARG001
    datadir: Path,
    tmp_path: Path,
) -> None:
    """Happy-path: extract AFs from local test variant and sample data."""
    in_sites = str(datadir / "chr1_100001_200000_0.01_seed42.ht")
    in_samples = str(datadir / "hgdp_1kg_sample_metadata.ht")
    out_va = str(tmp_path / "va.ht")
    out_sa = str(tmp_path / "sa.ht")

    with patch("divref.tools.extract_gnomad_afs.hl.init"):
        extract_gnomad_afs(
            in_gnomad_sites_table=in_sites,
            in_gnomad_hgdp_sample_data=in_samples,
            out_variant_annotation_table=out_va,
            out_sample_metadata=out_sa,
        )

    # --- variant annotation table ---
    va = hl.read_table(out_va)
    va_count = va.count()
    assert va_count > 0, "Expected at least one variant after filtering"

    # Globals should contain the requested populations
    pops = hl.eval(va.globals.pops)
    assert pops == DEFAULT_POPULATIONS

    # Each row should have pop_freqs with one entry per population
    first_row = va.head(1).collect()[0]
    assert len(first_row.pop_freqs) == len(DEFAULT_POPULATIONS)
    assert all(hasattr(f, "AF") for f in first_row.pop_freqs)

    # All retained variants must exceed freq_threshold in at least one population
    min_max_af = va.aggregate(hl.agg.min(hl.max(va.pop_freqs.map(lambda x: x.AF))))
    assert min_max_af > 0.001

    # --- sample metadata table ---
    sa = hl.read_table(out_sa)
    sa_count = sa.count()
    assert sa_count > 0, "Expected at least one sample"

    # Each sample should have a pop field
    first_sample = sa.head(1).collect()[0]
    assert hasattr(first_sample, "pop")
    assert isinstance(first_sample.pop, str)
