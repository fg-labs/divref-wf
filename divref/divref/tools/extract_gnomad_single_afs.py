"""Tool to extract gnomAD variant and sample frequency data for the DivRef pipeline."""

from enum import StrEnum
from enum import unique
from pathlib import Path

import hail as hl

from divref import defaults
from divref.hail import hail_init
from divref.haplotype import to_hashable_items


@unique
class GnomadVersion(StrEnum):
    """Mapping of gnomAD versions to GCS sites tables."""

    JOINT_41 = "gs://gcp-public-data--gnomad/release/4.1/ht/joint/gnomad.joint.v4.1.sites.ht"
    HGDP_1KG_312 = "gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_variant_annotations.ht"


def extract_gnomad_single_afs(
    *,
    gnomad_version: GnomadVersion,
    out_sites_hail_table: Path,
    out_sites_tsv: Path,
    contig: str,
    freq_threshold: float = 0,
    populations: list[str] = defaults.POPULATIONS,
    reference_genome: str = defaults.REFERENCE_GENOME,
    gcs_credentials_path: Path = Path("~/.config/gcloud/application_default_credentials.json"),
) -> None:
    """
    Extract gnomAD variant and sample frequency data for downstream pipeline tools.

    Reads a gnomAD sites table and filters to variants above the frequency threshold in at least one
    population. Writes two outputs: a Hail table at `out_sites_hail_table` for downstream pipeline
    tools, and a flat TSV at `out_sites_tsv` with columns `variant` (contig:pos:ref:alt) and one
    allele-frequency column per population.

    Args:
        out_sites_hail_table: Output path for the Hail table.
        out_sites_tsv: Output path for the TSV file.
        gnomad_version: gnomAD sites table to use - the schema varies.
        contig: Contig to extract sites.
        freq_threshold: Minimum allele frequency in any population to retain a variant. Defaults 0,
            all variants returned.
        populations: List of population codes to extract frequencies for.
        reference_genome: Reference genome to use. Defaults to "GRCh38".
        gcs_credentials_path: Path to GCS default credentials JSON file.
    """
    hail_init(gcs_credentials_path.expanduser())

    ht: hl.Table
    match gnomad_version:
        case GnomadVersion.JOINT_41:
            ht = extract_from_joint_41(
                contig=contig,
                freq_threshold=freq_threshold,
                populations=populations,
                reference_genome=reference_genome,
            )
        case GnomadVersion.HGDP_1KG_312:
            ht = extract_from_hgdp_1kg_312(
                contig=contig,
                freq_threshold=freq_threshold,
                populations=populations,
                reference_genome=reference_genome,
            )
        case _:
            raise ValueError(f"Unrecognized gnomAD sites table version: {gnomad_version}")

    ht.naive_coalesce(64).write(str(out_sites_hail_table), overwrite=True)

    ht.select(
        variant=hl.format(
            "%s:%s:%s:%s",
            ht.locus.contig,
            hl.str(ht.locus.position),
            ht.alleles[0],
            ht.alleles[1],
        ),
        **{pop: ht.pop_freqs[i].AF for i, pop in enumerate(populations)},
    ).export(str(out_sites_tsv))


def extract_from_joint_41(
    contig: str,
    freq_threshold: float,
    populations: list[str],
    reference_genome: str,
) -> hl.Table:
    """Use the gnomAD 4.1 joint sites schema."""
    va_all = hl.read_table(GnomadVersion.JOINT_41.value)
    interval = hl.parse_locus_interval(contig, reference_genome=reference_genome)
    va = hl.filter_intervals(va_all, [interval])

    freq_meta = va.globals.joint_globals.freq_meta.collect()[0]
    map_to_index = {to_hashable_items(x): i for i, x in enumerate(freq_meta)}

    pop_indices = []
    for pop in populations:
        idx = map_to_index.get(to_hashable_items({"group": "adj", "gen_anc": pop}))
        if idx is None:
            raise ValueError(f"Population {pop!r} not found in gnomAD frequency metadata")
        pop_indices.append(idx)

    va = va.select_globals(pops=populations)
    va = va.select(pop_freqs=hl.literal(pop_indices).map(lambda i: va.joint.freq[i]))
    va = va.filter(hl.any(lambda x: x.AF > freq_threshold, va.pop_freqs))
    va = va.key_by()

    return va.select("locus", "alleles", "pop_freqs")


def extract_from_hgdp_1kg_312(
    contig: str,
    freq_threshold: float,
    populations: list[str],
    reference_genome: str,
) -> hl.Table:
    """Use the gnomAD 3.1.2 HGDP+1KG sites schema."""
    va_all = hl.read_table(GnomadVersion.HGDP_1KG_312.value)
    interval = hl.parse_locus_interval(contig, reference_genome=reference_genome)
    va = hl.filter_intervals(va_all, [interval])

    freq_meta = va.globals.gnomad_freq_meta.collect()[0]
    map_to_index = {to_hashable_items(x): i for i, x in enumerate(freq_meta)}

    pop_indices = []
    for pop in populations:
        idx = map_to_index.get(to_hashable_items({"group": "adj", "pop": pop}))
        if idx is None:
            raise ValueError(f"Population {pop!r} not found in gnomAD frequency metadata")
        pop_indices.append(idx)

    va = va.select_globals(pops=populations)
    va = va.select(pop_freqs=hl.literal(pop_indices).map(lambda i: va.gnomad_freq[i]))
    va = va.filter(hl.any(lambda x: x.AF > freq_threshold, va.pop_freqs))
    va = va.key_by()

    return va.select("locus", "alleles", "pop_freqs")
