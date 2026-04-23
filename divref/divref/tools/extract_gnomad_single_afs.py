"""Tool to extract gnomAD variant and sample frequency data for the DivRef pipeline."""

import operator
from dataclasses import dataclass
from enum import StrEnum
from enum import unique
from pathlib import Path

import hail as hl
from fgpyo.io import assert_path_is_writable

from divref import defaults
from divref.hail import hail_init
from divref.haplotype import to_hashable_items


@unique
class GnomadVersion(StrEnum):
    """Mapping of gnomAD versions to GCS sites tables."""

    JOINT_41 = "gs://gcp-public-data--gnomad/release/4.1/ht/joint/gnomad.joint.v4.1.sites.ht"
    GENOMES_312 = (
        "gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.sites.ht"
    )
    HGDP_1KG_312 = "gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_variant_annotations.ht"


@dataclass(frozen=True)
class _GnomadSchema:
    """Per-version schema config for reading population AFs from a gnomAD sites table."""

    freq_meta_path: str  # dot-path into table globals, e.g. "joint_globals.freq_meta"
    row_freq_field: str  # dot-path on table row, e.g. "joint.freq"
    pop_key: str  # key in freq_meta entries for the population code, e.g. "gen_anc" or "pop"


_GNOMAD_SCHEMA: dict[GnomadVersion, _GnomadSchema] = {
    GnomadVersion.JOINT_41: _GnomadSchema(
        freq_meta_path="joint_globals.freq_meta",
        row_freq_field="joint.freq",
        pop_key="gen_anc",
    ),
    GnomadVersion.GENOMES_312: _GnomadSchema(
        freq_meta_path="freq_meta",
        row_freq_field="freq",
        pop_key="pop",
    ),
    GnomadVersion.HGDP_1KG_312: _GnomadSchema(
        freq_meta_path="gnomad_freq_meta",
        row_freq_field="gnomad_freq",
        pop_key="pop",
    ),
}


def extract_gnomad_single_afs(
    *,
    gnomad_version: GnomadVersion,
    contig: str,
    freq_threshold: float = 0,
    no_apply_filters: bool = False,
    populations: list[str] = defaults.POPULATIONS,
    reference_genome: str = defaults.REFERENCE_GENOME,
    out_sites_hail_table: Path | None = None,
    out_sites_tsv: Path | None = None,
    gcs_credentials_path: Path = Path("~/.config/gcloud/application_default_credentials.json"),
) -> None:
    """
    Extract gnomAD variant and sample frequency data for downstream pipeline tools.

    Reads a gnomAD sites table and filters to variants above the frequency threshold in at least one
    population. Writes up to two outputs: a Hail table at `out_sites_hail_table` for downstream
    pipeline tools, and a flat TSV at `out_sites_tsv` with columns `variant` (contig:pos:ref:alt)
    and one allele-frequency column per population.

    At least one of `out_sites_hail_table` or `out_sites_tsv` must be defined.

    Args:
        gnomad_version: gnomAD sites table to use - the schema varies per version.
        contig: Contig to extract sites.
        freq_threshold: Minimum allele frequency in any population to retain a variant. Defaults 0,
            all variants returned.
        no_apply_filters: If set, don't apply any of the variant filters (e.g. VQSR, AC0)
        populations: List of population codes to extract frequencies for.
        reference_genome: Reference genome to use. Defaults to "GRCh38".
        out_sites_hail_table: Output path for the Hail table. Optional.
        out_sites_tsv: Output path for the TSV file. Optional.
        gcs_credentials_path: Path to GCS default credentials JSON file.
    """
    if out_sites_hail_table is None and out_sites_tsv is None:
        raise ValueError("At least one of out_sites_hail_table or out_sites_tsv must be provided")

    # validate output paths before starting Hail
    if out_sites_tsv is not None:
        assert_path_is_writable(out_sites_tsv)

    hail_init(gcs_credentials_path.expanduser())

    schema = _GNOMAD_SCHEMA[gnomad_version]

    va_all = hl.read_table(gnomad_version.value)
    interval = hl.parse_locus_interval(contig, reference_genome=reference_genome)
    va = hl.filter_intervals(va_all, [interval])

    freq_meta = operator.attrgetter(schema.freq_meta_path)(va.globals).collect()[0]
    map_to_index = {to_hashable_items(x): i for i, x in enumerate(freq_meta)}

    pop_indices = []
    for pop in populations:
        idx = map_to_index.get(to_hashable_items({"group": "adj", schema.pop_key: pop}))
        if idx is None:
            raise ValueError(f"Population {pop!r} not found in gnomAD frequency metadata")
        pop_indices.append(idx)

    if not no_apply_filters:
        if gnomad_version is GnomadVersion.JOINT_41:
            # exome and genome filters are separate
            va_exome = va.filter(hl.coalesce(hl.len(va.exomes.filters) == 0, True))
            va_genome = va.filter(hl.coalesce(hl.len(va.genomes.filters) == 0, True))
            va = va_exome.union(va_genome)
        elif gnomad_version in [GnomadVersion.GENOMES_312, GnomadVersion.HGDP_1KG_312]:
            # Some filter sets are {} and some are NA; treat NA as passing.
            va = va.filter(hl.coalesce(hl.len(va.filters) == 0, True))

    va = va.select_globals(pops=populations)
    row_freq = operator.attrgetter(schema.row_freq_field)(va)
    va = va.select(pop_freqs=hl.literal(pop_indices).map(lambda i: row_freq[i]))
    if freq_threshold > 0:
        va = va.filter(hl.any(lambda x: x.AF > freq_threshold, va.pop_freqs))
    va = va.key_by()
    va = va.select("locus", "alleles", "pop_freqs")

    if out_sites_hail_table is not None:
        va.naive_coalesce(64).write(str(out_sites_hail_table), overwrite=True)

    if out_sites_tsv is not None:
        va.select(
            variant=hl.format(
                "%s:%s:%s:%s",
                va.locus.contig,
                hl.str(va.locus.position),
                va.alleles[0],
                va.alleles[1],
            ),
            **{pop: va.pop_freqs[i].AF for i, pop in enumerate(populations)},
        ).export(str(out_sites_tsv))
