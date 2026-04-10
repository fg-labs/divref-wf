"""Tool to extract gnomAD variant and sample frequency data for the DivRef pipeline."""

import hail as hl

from divref.haplotype import HailPath
from divref.haplotype import to_hashable_items

_GNOMAD_SITES_TABLE_PATH = (
    "gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/"
    "gnomad.genomes.v3.1.2.hgdp_1kg_subset_variant_annotations.ht"
)
_GNOMAD_HGDP_SAMPLE_DATA_PATH = (
    "gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/"
    "gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.ht"
)
_HGDP_REPRESENTED_POPS = ["afr", "amr", "eas", "sas", "nfe"]


def extract_gnomad_afs(
    *,
    output_file_va: HailPath,
    output_file_sa: HailPath,
    freq_threshold: float = 0.001,
    gnomad_sites_table_path: HailPath = _GNOMAD_SITES_TABLE_PATH,
    gnomad_hgdp_sample_data_path: HailPath = _GNOMAD_HGDP_SAMPLE_DATA_PATH,
    hgdp_represented_pops: list[str] = _HGDP_REPRESENTED_POPS,
) -> None:
    """
    Extract gnomAD variant and sample frequency data for downstream pipeline tools.

    Reads the gnomAD v3.1.2 HGDP/1KG subset, extracts per-population allele frequencies
    for the specified populations, filters to variants above the frequency threshold in at
    least one population, and writes compact variant annotation and sample metadata tables.

    Args:
        output_file_va: Output path for the variant annotation Hail table.
        output_file_sa: Output path for the sample metadata Hail table.
        freq_threshold: Minimum allele frequency in any population to retain a variant.
        gnomad_sites_table_path: Path to the gnomAD v3.1.2 HGDP/1KG variant annotations table.
        gnomad_hgdp_sample_data_path: Path to the gnomAD v3.1.2 HGDP/1KG sample metadata table.
        hgdp_represented_pops: List of population codes to extract frequencies for.
    """
    hl.init()

    va = hl.read_table(gnomad_sites_table_path)

    freq_meta = va.globals.gnomad_freq_meta[:15].collect()[0]
    map_to_index = {to_hashable_items(x): i for i, x in enumerate(freq_meta)}

    pop_indices = []
    for pop in hgdp_represented_pops:
        idx = map_to_index.get(to_hashable_items({"group": "adj", "pop": pop}))
        if idx is None:
            raise ValueError(f"Population {pop!r} not found in gnomAD frequency metadata")
        pop_indices.append(idx)

    # Some filter sets are {} and some are NA; treat NA as passing.
    va = va.filter(hl.coalesce(hl.len(va.filters) == 0, True))

    va = va.select_globals(pops=hgdp_represented_pops)
    va = va.select(pop_freqs=hl.literal(pop_indices).map(lambda i: va.gnomad_freq[i]))
    va = va.filter(hl.any(lambda x: x.AF > freq_threshold, va.pop_freqs))
    va.naive_coalesce(64).write(output_file_va, overwrite=True)

    sa = hl.read_table(gnomad_hgdp_sample_data_path).select_globals()
    sa = sa.select(pop=sa.gnomad_population_inference.pop)
    sa.naive_coalesce(4).write(output_file_sa, overwrite=True)
