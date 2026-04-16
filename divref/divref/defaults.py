from typing import Final

from divref.alias import HailPath

GNOMAD_HGDP_1KG_VARIANT_ANNOTATION_HAIL_TABLE: Final[HailPath] = (
    "gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_variant_annotations.ht"  # noqa: E501
)
"""HGDP+1KG individual level genotypes."""

GNOMAD_HGDP_1KG_SAMPLE_METADATA_HAIL_TABLE: Final[HailPath] = (
    "gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.ht"  # noqa: E501
)
"""HGDP+1KG sample metadata."""

POPULATIONS: list[str] = ["afr", "amr", "eas", "sas", "nfe"]
"""Default HGDP+1KG populations."""
