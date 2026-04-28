"""Tool to build DivRef DuckDB index from haplotype Hail tables."""

import json
import logging
import os
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
) -> None:
    """
    Convert per-chr haplotype and gnomAD variant Hail tables into a searchable DuckDB index.

    Reads and annotates the haplotype and gnomAD variant tables one contig at a time, sorting by
    position within each contig.

    Merges all the per-contig tables together, assigns sequence IDs, generates sequence strings with
    flanking reference context.

    Writes a DuckDB index file.

    Args:
        in_table_pairs_tsv: Path to a TSV file with fields 'contig', 'haplotype_table_path', and
            'sites_table_path'.
        reference_fasta: Path to the indexed reference FASTA for sequence extraction.
        window_size: Window size used when generating haplotypes; used as the context size when
            constructing sequence strings and stored in the index.
        output_base: Base path for output files. Writes {output_base}.haplotypes.tsv.bgz and
            {output_base}.haplotypes.index.duckdb.
        version: Version identifier embedded in sequence IDs (e.g. "1.0").
        reference_genome: Reference genome to use. Defaults to "GRCh38".
        tmp_dir: Temporary directory for Hail checkpoint files.
    """
    assert_path_is_readable(in_table_pairs_tsv)
    assert_path_is_readable(reference_fasta)
    assert_path_is_readable(reference_fasta.with_suffix(".fai"))
    assert_directory_exists(tmp_dir)

    # fail fast if output files are not writeable
    out_table_file: Path = Path(f"{str(output_base)}.haplotypes_gnomad_merge.tsv.bgz")
    out_duckdb_file: Path = Path(f"{str(output_base)}.haplotypes_gnomad_merge.index.duckdb")
    assert_path_is_writable(out_table_file)
    assert_path_is_writable(out_duckdb_file)

    table_pairs: list[TablePair] = list(TablePair.read(in_table_pairs_tsv))

    # fail fast on input Hail tables
    for table_pair in table_pairs:
        assert_directory_exists(table_pair.haplotype_table_path)
        assert_directory_exists(table_pair.sites_table_path)

    hl.init(tmp_dir=str(tmp_dir))

    pops_legend: list[str] = hl.read_table(str(table_pairs[0].sites_table_path)).pops.collect()[0]

    ht = build_sequences_table(
        table_pairs=table_pairs,
        reference_fasta=reference_fasta,
        reference_genome=reference_genome,
        window_size=window_size,
        version=version,
    )
    df = export_sequences_table_to_dataframe(
        ht=ht, out_file=out_table_file, pops_legend=pops_legend
    )
    write_duckdb_index(
        df=df,
        out_file=out_duckdb_file,
        window_size=window_size,
        pops_legend=pops_legend,
        version=version,
    )


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
    argmax_pop = hl.argmax(va.gnomad_freqs.map(lambda x: hl.max(x.AF)))
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


def build_sequences_table(
    table_pairs: list[TablePair],
    reference_fasta: Path,
    reference_genome: str,
    window_size: int,
    version: str,
) -> hl.Table:
    """
    Build the "sequences" table input with sequences and variant strings.

    Prepares and annotates the haplotype and gnomAD tables, assigns sequence IDs, and generates
    sequences with flanking reference context.

    Args:
        table_pairs: List of table pairs by contig. Sequences will be sorted by reference genome
            start position within each contig; contigs will be ordered as input.
        reference_fasta: Path to the reference FASTA.
        reference_genome: Reference genome to use.
        window_size: Context size for sequence construction and haplotype splitting.
        version: Version identifier for sequence IDs.

    Returns:
        Hail table with sequences, coordinates, and variant strings annotated.
    """
    hl.get_reference(reference_genome).add_sequence(str(reference_fasta))

    sequences_tables: list[hl.Table] = []
    for table_pair in table_pairs:
        hgdp_haplotypes_ht: hl.Table = build_hgdp_haplotype_table_entries(
            haplotypes_table_path=table_pair.haplotype_table_path, window_size=window_size
        )
        gnomad_variants_ht: hl.Table = build_gnomad_variant_table_entries(
            sites_table_path=table_pair.sites_table_path
        )
        contig_seq_ht: hl.Table = hgdp_haplotypes_ht.union(gnomad_variants_ht, unify=True)
        sequences_tables.append(contig_seq_ht)

    seq_ht: hl.Table = sequences_tables[0].union(*sequences_tables[1:])

    seq_ht = seq_ht.rename({
        "max_empirical_AF": "popmax_empirical_AF",
        "max_empirical_AC": "popmax_empirical_AC",
    })

    seq_ht = seq_ht.annotate(
        min_locus=hl.sorted(seq_ht.variants, key=lambda v: v.locus.position)[0].locus
    )
    seq_ht = seq_ht.order_by(seq_ht.min_locus).drop("min_locus")
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
        sequence_id=hl.str(f"DR-{version}-") + hl.str(seq_ht.idx),
        n_variants=hl.len(seq_ht.variants),
    ).drop("idx")

    return seq_ht


def export_sequences_table_to_dataframe(
    ht: hl.Table,
    out_file: Path,
    pops_legend: list[str],
) -> polars.DataFrame:
    """
    Export the sequences Hail table to a TSV file and return it as a polars DataFrame.

    Args:
        ht: Annotated haplotype/variant table with sequences and variant strings.
        out_file: Path for the output TSV file.
        pops_legend: Ordered list of population codes for frequency columns.

    Returns:
        Polars DataFrame read back from the exported TSV.
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

    schema_overrides: dict[str, type[polars.DataType]] = {
        "sequence_id": polars.String,
        **{f"gnomAD_AF_{pop}": polars.String for pop in pops_legend},
    }
    return polars.read_csv(
        out_file,
        separator="\t",
        schema_overrides=schema_overrides,
        null_values="null",
    )


def write_duckdb_index(
    df: polars.DataFrame,  # noqa: ARG001 — accessed by DuckDB via SQL `FROM df`
    out_file: Path,
    window_size: int,
    pops_legend: list[str],
    version: str,
) -> None:
    """
    Create a DuckDB index database from the haplotype DataFrame.

    Args:
        df: DataFrame with sequence and metadata columns.
        out_file: Path for the output .duckdb file.
        window_size: Window size to store in the index.
        pops_legend: Population legend list to store in the index.
        version: Version string to store in the index.
    """
    if os.path.exists(out_file):
        os.remove(out_file)
    con = duckdb.connect(str(out_file))
    con.execute("CREATE TABLE sequences AS SELECT * FROM df")
    con.execute("CREATE INDEX idx_sequence_id ON sequences(sequence_id)")
    con.execute("CREATE TABLE window_size AS SELECT ? AS window_size", [window_size])
    con.execute("CREATE TABLE pops_legend AS SELECT ? AS pops_legend", [json.dumps(pops_legend)])
    con.execute("CREATE TABLE VERSION AS SELECT ? AS version", [version])
    con.close()
