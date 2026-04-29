"""Tool to build DivRef DuckDB index from haplotype Hail tables."""

import json
import logging
from collections.abc import Iterator
from pathlib import Path

import duckdb
import hail as hl
import polars
from fgmetric import Metric
from fgpyo.io import assert_directory_exists
from fgpyo.io import assert_path_is_readable
from fgpyo.io import assert_path_is_writable

from divref import defaults
from divref.haplotype import get_haplo_sequence
from divref.haplotype import haplo_coordinates
from divref.haplotype import split_haplotypes

logger = logging.getLogger(__name__)


class TablePair(Metric):
    """
    Helper class to link a pair of tables for the same contig.

    Attributes:
        contig: Contig name.
        haplotype_table: HGDP haplotypes Hail table.
        sites_table: gnomAD variant Hail table.
    """

    contig: str
    haplotype_table_path: Path
    sites_table_path: Path


def create_duckdb_index(
    *,
    in_table_pairs_tsv: Path,
    reference_fasta: Path,
    window_size: int,
    output_base: Path,
    version: str,
    reference_genome: str = defaults.REFERENCE_GENOME,
    tmp_dir: Path = Path("/tmp"),
    polars_chunk_size: int = 100_000,
    retain_per_contig_tsvs: bool = False,
    force: bool = False,
) -> None:
    """
    Convert per-chr haplotype and gnomAD variant Hail tables into a searchable DuckDB index.

    Streams sequences one contig at a time, then sub-divides each contig at the polars read step:
    each per-contig TSV is read back in batches of `polars_chunk_size` rows, and each batch is
    appended to the DuckDB `sequences` table. The first batch creates the table; subsequent
    batches `INSERT INTO` it. Sequence IDs are assigned with a running offset so they remain
    unique across contigs and batches.

    Args:
        in_table_pairs_tsv: Path to a TSV file with fields 'contig', 'haplotype_table_path', and
            'sites_table_path'.
        reference_fasta: Path to the indexed reference FASTA for sequence extraction.
        window_size: Window size used when generating haplotypes; used as the context size when
            constructing sequence strings and stored in the index.
        output_base: Base path for output. Writes
            `{output_base}.haplotypes_gnomad_merge.index.duckdb` and, when
            `retain_per_contig_tsvs` is True, one
            `{output_base}.haplotypes_gnomad_merge.{contig}.tsv.bgz` per contig.
        version: Version identifier embedded in sequence IDs (e.g. "1.0").
        reference_genome: Reference genome to use. Defaults to "GRCh38".
        tmp_dir: Temporary directory for Hail checkpoints and (when not retained) per-contig
            intermediate TSVs.
        polars_chunk_size: Maximum number of rows per polars read batch. Bounds the in-process
            DataFrame size when streaming each per-contig TSV into DuckDB.
        retain_per_contig_tsvs: If True, write per-contig TSVs alongside the duckdb output rather
            than into `tmp_dir`.
        force: If True, overwrite an existing duckdb output. Otherwise raise FileExistsError.
    """
    assert_path_is_readable(in_table_pairs_tsv)
    assert_path_is_readable(reference_fasta)
    assert_path_is_readable(reference_fasta.with_suffix(".fai"))
    assert_directory_exists(tmp_dir)

    out_duckdb_file: Path = Path(f"{str(output_base)}.haplotypes_gnomad_merge.index.duckdb")
    if out_duckdb_file.exists():
        if not force:
            raise FileExistsError(
                f"DuckDB output already exists at {out_duckdb_file}. Pass --force to overwrite."
            )
        out_duckdb_file.unlink()
    assert_path_is_writable(out_duckdb_file)

    table_pairs: list[TablePair] = list(TablePair.read(in_table_pairs_tsv))

    # fail fast on input Hail tables
    for table_pair in table_pairs:
        assert_directory_exists(table_pair.haplotype_table_path)
        assert_directory_exists(table_pair.sites_table_path)

    # determine per-contig TSV paths and fail fast on writability when retained
    per_contig_tsv_dir: Path = output_base.parent if retain_per_contig_tsvs else tmp_dir
    per_contig_tsvs: dict[str, Path] = {
        tp.contig: per_contig_tsv_dir
        / f"{output_base.name}.haplotypes_gnomad_merge.{tp.contig}.tsv.bgz"
        for tp in table_pairs
    }
    if retain_per_contig_tsvs:
        for tsv_path in per_contig_tsvs.values():
            assert_path_is_writable(tsv_path)

    hl.init(tmp_dir=str(tmp_dir))

    pops_legend: list[str] = hl.read_table(str(table_pairs[0].sites_table_path)).pops.collect()[0]
    hl.get_reference(reference_genome).add_sequence(str(reference_fasta))

    with duckdb.connect(str(out_duckdb_file)) as conn:
        sequence_id_offset: int = 0
        created_table: bool = False
        for table_pair in table_pairs:
            contig_seq_ht = build_contig_sequences_table(
                table_pair=table_pair,
                window_size=window_size,
                version=version,
                sequence_id_offset=sequence_id_offset,
            )
            contig_tsv: Path = per_contig_tsvs[table_pair.contig]
            export_sequences_table_to_tsv(
                ht=contig_seq_ht, out_file=contig_tsv, pops_legend=pops_legend
            )

            contig_rows: int = 0
            for df in iter_dataframe_chunks(
                tsv=contig_tsv, pops_legend=pops_legend, chunk_size=polars_chunk_size
            ):
                if not created_table:
                    conn.execute("CREATE TABLE sequences AS SELECT * FROM df")
                    created_table = True
                else:
                    conn.execute("INSERT INTO sequences SELECT * FROM df")
                contig_rows += df.height
                sequence_id_offset += df.height

            if not retain_per_contig_tsvs and contig_tsv.exists():
                contig_tsv.unlink()

            logger.info(
                f"Appended {contig_rows} rows for contig {table_pair.contig} "
                f"(running total: {sequence_id_offset})"
            )

        conn.execute("CREATE INDEX idx_sequence_id ON sequences(sequence_id)")
        conn.execute("CREATE TABLE window_size AS SELECT ? AS window_size", [window_size])
        conn.execute(
            "CREATE TABLE pops_legend AS SELECT ? AS pops_legend", [json.dumps(pops_legend)]
        )
        conn.execute("CREATE TABLE VERSION AS SELECT ? AS version", [version])


def build_hgdp_haplotype_table_entries(
    haplotypes_table_path: Path,
    window_size: int,
) -> hl.Table:
    """
    Build HGDP_haplotype entries for the "sequences" table.

    Reads the haplotype table, splits the haplotypes by window size, and annotates with source and
    population frequencies.

    Args:
        haplotypes_table_path: Path to the computed haplotypes Hail table.
        window_size: Context size for sequence construction and haplotype splitting.

    Returns:
        Hail table with added sequences and variant strings.
    """
    # Read the table and remove keys
    ht = hl.read_table(str(haplotypes_table_path)).key_by()
    count_orig: int = ht.count()
    logger.info(f"Haplotype table {haplotypes_table_path} contains {count_orig} unique haplotypes.")

    # Split haplotypes by window size
    ht = split_haplotypes(ht, window_size)
    ht = ht.key_by("haplotype").distinct().key_by().drop("haplotype")
    count_after_splitting: int = ht.count()
    logger.info(
        f"{count_after_splitting} unique haplotypes remaining after splitting at "
        f"window size {window_size}"
    )

    # Annotate
    ht = ht.annotate(
        source="HGDP_haplotype",
        all_pop_freqs=ht.all_pop_freqs.map(
            lambda x: hl.struct(pop=x.pop, empirical_AC=x.empirical_AC, empirical_AF=x.empirical_AF)
        ),
    )

    return ht


def build_gnomad_variant_table_entries(sites_table_path: Path) -> hl.Table:
    """
    Build gnomAD_variant entries for the "sequences" table.

    Reads the gnomAD table and annotates entries to match the HGDP_haplotype entries.

    Args:
        sites_table_path: Path to the gnomAD variant annotations Hail table.

    Returns:
        Tuple of (checkpointed Hail table, population legend list).
    """
    va = hl.read_table(str(sites_table_path))
    count_orig: int = va.count()
    logger.info(f"Variant table {sites_table_path} contains {count_orig} variants.")

    va = va.rename({"pop_freqs": "gnomad_freqs"})
    va = va.key_by()
    argmax_pop = hl.argmax(va.gnomad_freqs.map(lambda x: x.AF))
    va = va.select(
        max_pop=argmax_pop,
        max_empirical_AF=va.gnomad_freqs[argmax_pop].AF,
        fraction_phased=1.0,
        estimated_gnomad_AF=va.gnomad_freqs[argmax_pop].AF,
        max_empirical_AC=va.gnomad_freqs[argmax_pop].AC,
        all_pop_freqs=hl.range(hl.len(va.gnomad_freqs)).map(
            lambda i: hl.struct(
                pop=i,
                empirical_AC=va.gnomad_freqs[i].AC,
                empirical_AF=va.gnomad_freqs[i].AF,
            )
        ),
        source="gnomAD_variant",
        variants=[hl.struct(locus=va.locus, alleles=va.alleles)],
        gnomad_freqs=[va.gnomad_freqs],
    )
    return va


def build_contig_sequences_table(
    *,
    table_pair: TablePair,
    window_size: int,
    version: str,
    sequence_id_offset: int,
) -> hl.Table:
    """
    Build the per-contig sequences hail table with sequences, coordinates, and IDs.

    Reads the HGDP haplotype + gnomAD sites tables for one contig, unions them, sorts by genomic
    position, and applies the same per-row annotations as the cross-contig table. Sequence IDs are
    offset by `sequence_id_offset` so they remain unique across contigs.

    Args:
        table_pair: Per-contig pair of haplotype + gnomAD sites table paths.
        window_size: Context size for sequence construction and haplotype splitting.
        version: Version identifier for sequence IDs.
        sequence_id_offset: Number of rows already written for prior contigs; added to this
            contig's local index to produce a globally unique sequence ID.

    Returns:
        Hail table with sequences, coordinates, and variant strings annotated.
    """
    hgdp_haplotypes_ht: hl.Table = build_hgdp_haplotype_table_entries(
        haplotypes_table_path=table_pair.haplotype_table_path, window_size=window_size
    )
    gnomad_variants_ht: hl.Table = build_gnomad_variant_table_entries(
        sites_table_path=table_pair.sites_table_path
    )
    seq_ht: hl.Table = hgdp_haplotypes_ht.union(gnomad_variants_ht, unify=True)

    seq_ht = seq_ht.rename({
        "max_empirical_AF": "popmax_empirical_AF",
        "max_empirical_AC": "popmax_empirical_AC",
    })

    seq_ht = seq_ht.annotate(
        min_pos=hl.sorted(seq_ht.variants, key=lambda v: v.locus.position)[0].locus.position
    )
    seq_ht = seq_ht.order_by(seq_ht.min_pos).drop("min_pos")
    seq_ht = seq_ht.add_index()
    coords = haplo_coordinates(window_size, seq_ht.variants)
    seq_ht = seq_ht.annotate(
        sequence=get_haplo_sequence(window_size, seq_ht.variants),
        contig=seq_ht.variants[0].locus.contig,
        start=coords.start,
        end=coords.end,
    )
    seq_ht = seq_ht.annotate(variant_strs=seq_ht.variants.map(lambda x: hl.variant_str(x)))
    seq_ht = seq_ht.annotate(
        sequence_length=hl.len(seq_ht.sequence),
        sequence_id=hl.str(f"DR-{version}-") + hl.str(seq_ht.idx + sequence_id_offset),
        n_variants=hl.len(seq_ht.variants),
    ).drop("idx")

    return seq_ht


def export_sequences_table_to_tsv(
    ht: hl.Table,
    out_file: Path,
    pops_legend: list[str],
) -> None:
    """
    Export the sequences Hail table to a single bgz-compressed TSV.

    Args:
        ht: Annotated haplotype/variant table with sequences and variant strings.
        out_file: Path for the output TSV file.
        pops_legend: Ordered list of population codes for frequency columns.
    """
    ht.select(
        "sequence",
        "sequence_length",
        "sequence_id",
        "n_variants",
        "contig",
        "start",
        "end",
        "popmax_empirical_AF",
        "popmax_empirical_AC",
        "estimated_gnomad_AF",
        "fraction_phased",
        "source",
        max_pop=hl.literal(pops_legend)[ht.max_pop],
        variants=hl.delimit(ht.variant_strs, ","),
        **{
            f"gnomAD_AF_{pop}": hl.delimit(
                ht.gnomad_freqs.map(lambda x, _i=i: hl.format("%.5f", x[_i].AF)), ","
            )
            for i, pop in enumerate(pops_legend)
        },
    ).export(str(out_file))


def iter_dataframe_chunks(
    *,
    tsv: Path,
    pops_legend: list[str],
    chunk_size: int,
) -> Iterator[polars.DataFrame]:
    """
    Yield polars DataFrames of up to `chunk_size` rows from a sequences TSV.

    The `sequence_id` and `gnomAD_AF_*` columns are explicitly typed as strings so that
    schema inference cannot misread comma-delimited per-variant AFs as floats.

    Args:
        tsv: Path to the sequences TSV (bgz-compressed).
        pops_legend: Ordered list of population codes used to name `gnomAD_AF_{pop}` columns.
        chunk_size: Maximum rows per yielded DataFrame.

    Yields:
        Polars DataFrame batches read from `tsv`.
    """
    schema_overrides: dict[str, type[polars.DataType]] = {
        "sequence_id": polars.String,
        **{f"gnomAD_AF_{pop}": polars.String for pop in pops_legend},
    }
    lf = polars.scan_csv(
        tsv,
        separator="\t",
        schema_overrides=schema_overrides,
        null_values="null",
    )
    for df in lf.collect_batches(chunk_size=chunk_size):
        if df.height > 0:
            yield df
