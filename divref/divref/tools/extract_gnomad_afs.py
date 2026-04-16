"""Tool to extract gnomAD variant and sample frequency data for the DivRef pipeline."""

from pathlib import Path

import hail as hl

from divref import defaults
from divref.alias import HailPath
from divref.hail import hail_init
from divref.haplotype import to_hashable_items


def extract_gnomad_afs(
    *,
    in_gnomad_sites_table: HailPath,
    out_variant_annotation_table: HailPath,
    freq_threshold: float = 0.001,
    populations: list[str] = defaults.POPULATIONS,
    gcs_credentials_path: Path = Path("~/.config/gcloud/application_default_credentials.json"),
) -> None:
    """
    Extract gnomAD variant and sample frequency data for downstream pipeline tools.

    Reads the gnomAD v3.1.2 HGDP/1KG subset, extracts per-population allele frequencies
    for the specified populations, filters to variants above the frequency threshold in at
    least one population, and writes compact variant annotation and sample metadata tables.

    Args:
        in_gnomad_sites_table: Path to the gnomAD HGDP/1KG sites table.
        out_variant_annotation_table: Output path for the variant annotation Hail table.
        freq_threshold: Minimum allele frequency in any population to retain a variant.
        populations: List of population codes to extract frequencies for.
        gcs_credentials_path: Path to GCS default credentials JSON file.
    """
    hail_init(gcs_credentials_path.expanduser())

    va = hl.read_table(in_gnomad_sites_table)

    freq_meta = va.globals.gnomad_freq_meta.collect()[0]
    map_to_index = {to_hashable_items(x): i for i, x in enumerate(freq_meta)}

    pop_indices = []
    for pop in populations:
        idx = map_to_index.get(to_hashable_items({"group": "adj", "pop": pop}))
        if idx is None:
            raise ValueError(f"Population {pop!r} not found in gnomAD frequency metadata")
        pop_indices.append(idx)

    # Some filter sets are {} and some are NA; treat NA as passing.
    va = va.filter(hl.coalesce(hl.len(va.filters) == 0, True))

    va = va.select_globals(pops=populations)
    va = va.select(pop_freqs=hl.literal(pop_indices).map(lambda i: va.gnomad_freq[i]))
    va = va.filter(hl.any(lambda x: x.AF > freq_threshold, va.pop_freqs))
    va.naive_coalesce(64).write(out_variant_annotation_table, overwrite=True)
