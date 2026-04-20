"""Tool to build DivRef FASTA sequences and DuckDB index from haplotype Hail tables."""

import json
import logging
import os
from pathlib import Path

import duckdb
import hail as hl
import polars
from fgpyo.io import assert_directory_exists
from fgpyo.io import assert_fasta_indexed
from fgpyo.io import assert_path_is_readable
from fgpyo.io import assert_path_is_writable

from divref.haplotype import get_haplo_sequence
from divref.haplotype import split_haplotypes

logger = logging.getLogger(__name__)


def build_haplotype_table(
    haplotypes_table_path: Path,
    gnomad_va_file: Path,
    reference_fasta: Path,
    window_size: int,
    frequency_cutoff: float,
    merge: bool,
    version_str: str,
    tmp_dir: Path,
) -> tuple[hl.Table, list[str]]:
    """
    Build the annotated haplotype table with sequences and variant strings.

    Reads the haplotype and gnomAD tables, filters by estimated frequency, optionally
    merges single-variant gnomAD sites, assigns sequence IDs, and generates haplotype
    sequences with flanking reference context.

    Args:
        haplotypes_table_path: Path to the computed haplotypes Hail table.
        gnomad_va_file: Path to the gnomAD variant annotations Hail table.
        reference_fasta: Path to the GRCh38 reference FASTA.
        window_size: Context size for sequence construction and haplotype splitting.
        frequency_cutoff: Minimum estimated gnomAD AF for inclusion.
        merge: If True, include gnomAD single-variant sites above frequency_cutoff.
        version_str: Version identifier for sequence IDs.
        tmp_dir: Directory for Hail checkpoint files.

    Returns:
        Tuple of (checkpointed Hail table, population legend list).
    """
    ht = hl.read_table(str(haplotypes_table_path)).key_by()
    va = hl.read_table(str(gnomad_va_file))
    pops_legend: list[str] = va.pops.collect()[0]

    hl.get_reference("GRCh38").add_sequence(reference_fasta)

    logger.info(
        "Haplotype table contains %d unique haplotypes above frequency threshold", ht.count()
    )

    count_before = ht.count()
    ht = ht.filter(ht.min_variant_frequency > 0)
    count_after = ht.count()
    if count_after < count_before:
        logger.warning(
            "Removed %d haplotypes with min_variant_frequency <= 0", count_before - count_after
        )

    fraction_phased = ht.max_empirical_AF / ht.min_variant_frequency
    ht = ht.annotate(
        fraction_phased=fraction_phased,
        estimated_gnomad_AF=hl.min(
            ht.gnomad_freqs.map(lambda x: x[ht.max_pop].AF * fraction_phased)
        ),
    )
    ht = ht.filter(ht.estimated_gnomad_AF >= frequency_cutoff)
    ht = ht.rename({"max_empirical_AN": "max_empirical_AC"})
    ht = split_haplotypes(ht, window_size)
    ht = ht.key_by("haplotype").distinct().key_by().drop("haplotype")
    ht = ht.annotate(
        source="HGDP_haplotype",
        all_pop_freqs=ht.all_pop_freqs.map(
            lambda x: hl.struct(pop=x.pop, empirical_AC=x.empirical_AN, empirical_AF=x.empirical_AF)
        ),
    )
    logger.info(
        "After splitting at window size %d: %d unique haplotypes above frequency threshold",
        window_size,
        ht.count(),
    )

    if merge:
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
        va = va.filter(va.max_empirical_AF >= frequency_cutoff)
        ht = ht.union(va, unify=True)

    ht = ht.rename({
        "max_empirical_AF": "popmax_empirical_AF",
        "max_empirical_AC": "popmax_empirical_AC",
    })

    ht = ht.add_index()
    ht = ht.annotate(sequence=get_haplo_sequence(window_size, ht.variants))
    ht = ht.annotate(variant_strs=ht.variants.map(lambda x: hl.variant_str(x)))
    ht = ht.annotate(
        sequence_length=hl.len(ht.sequence),
        sequence_id=hl.str(f"DR-{version_str}-") + hl.str(ht.idx),
        n_variants=hl.len(ht.variants),
    ).drop("idx")

    file_suffix = ".haplotypes" if not merge else ".haplotypes_gnomad_merge"
    ht = ht.checkpoint(os.path.join(tmp_dir, f"{file_suffix}.ht"), overwrite=True)

    return ht, pops_legend


def export_ht_to_dataframe(
    ht: hl.Table,
    output_base: Path,
    file_suffix: str,
    pops_legend: list[str],
) -> polars.DataFrame:
    """
    Export the haplotype Hail table to a TSV file and return it as a polars DataFrame.

    Args:
        ht: Annotated haplotype table with sequences and variant strings.
        output_base: Base path for the output TSV file.
        file_suffix: Suffix to append before .tsv.bgz (e.g. ".haplotypes").
        pops_legend: Ordered list of population codes for frequency columns.

    Returns:
        Polars DataFrame read back from the exported TSV.
    """
    out_file_path: str = f"{str(output_base)}{file_suffix}.tsv.bgz"
    ht.select(
        "sequence",
        "sequence_length",
        "sequence_id",
        "n_variants",
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
    ).export(out_file_path)

    return polars.read_csv(
        out_file_path,
        separator="\t",
        schema_overrides={"sequence_id": polars.String},
    )


def write_fasta_files(
    df: polars.DataFrame,
    output_base: Path,
    file_suffix: str,
    split_contigs: bool,
) -> None:
    """
    Write FASTA file(s) from the haplotype DataFrame.

    Args:
        df: DataFrame with sequence and sequence_id columns.
        output_base: Base path for the output FASTA file(s).
        file_suffix: Suffix to append before .fasta (e.g. ".haplotypes").
        split_contigs: If True, write one FASTA file per chromosome.
    """
    if split_contigs:
        df = df.with_columns(contig=df["variants"].str.split(":").list.get(0))
        for chrom in df["contig"].unique().to_list():
            logger.info("Creating FASTA for chromosome %s", chrom)
            df2 = df.filter(df["contig"] == chrom)
            with open(f"{str(output_base)}{file_suffix}.{chrom}.fasta", "w") as fasta_out:
                for sequence, sequence_id in df2.select("sequence", "sequence_id").iter_rows():
                    fasta_out.write(f">{sequence_id}\n{sequence}\n")
    else:
        logger.info("Creating FASTA")
        with open(f"{str(output_base)}{file_suffix}.fasta", "w") as fasta_out:
            for sequence, sequence_id in df.select("sequence", "sequence_id").iter_rows():
                fasta_out.write(f">{sequence_id}\n{sequence}\n")


def create_duckdb_index(
    df: polars.DataFrame,  # noqa: ARG001 — accessed by DuckDB via SQL `FROM df`
    output_base: Path,
    file_suffix: str,
    window_size: int,
    pops_legend: list[str],
    version_str: str,
) -> None:
    """
    Create a DuckDB index database from the haplotype DataFrame.

    Args:
        df: DataFrame with sequence and metadata columns.
        output_base: Base path for the output .duckdb file.
        file_suffix: Suffix to append before .index.duckdb.
        window_size: Window size to store in the index.
        pops_legend: Population legend list to store in the index.
        version_str: Version string to store in the index.
    """
    duckdb_file = f"{str(output_base)}{file_suffix}.index.duckdb"
    if os.path.exists(duckdb_file):
        os.remove(duckdb_file)
    con = duckdb.connect(duckdb_file)
    con.execute("CREATE TABLE sequences AS SELECT * FROM df")
    con.execute("CREATE INDEX idx_sequence_id ON sequences(sequence_id)")
    con.execute("CREATE TABLE window_size AS SELECT ? AS window_size", [window_size])
    con.execute("CREATE TABLE pops_legend AS SELECT ? AS pops_legend", [json.dumps(pops_legend)])
    con.execute("CREATE TABLE VERSION AS SELECT ? AS version", [version_str])
    con.close()


def create_fasta_and_index(
    *,
    haplotypes_table_path: Path,
    gnomad_va_file: Path,
    reference_fasta: Path,
    window_size: int,
    output_base: Path,
    version_str: str,
    merge: bool = False,
    frequency_cutoff: float = 0.005,
    split_contigs: bool = False,
    tmp_dir: Path = Path("/tmp"),
) -> None:
    """
    Convert a haplotype Hail table into FASTA sequences and a searchable DuckDB index.

    Reads the haplotype table, filters by estimated gnomAD allele frequency, optionally
    merges in single gnomAD variants, assigns sequence IDs, generates sequence strings
    with flanking reference context, and writes FASTA and DuckDB index files for use by
    remap_divref.

    Args:
        haplotypes_table_path: Path to the Hail table of computed haplotypes
            (from compute_haplotypes).
        gnomad_va_file: Path to the gnomAD variant annotations Hail table
            (from extract_gnomad_afs).
        reference_fasta: Path to the GRCh38 reference FASTA for sequence extraction.
        window_size: Window size used when generating haplotypes; used as the context
            size when constructing sequence strings and stored in the index.
        output_base: Base path for output files. Writes {output_base}.haplotypes.tsv.bgz,
            {output_base}.haplotypes.fasta (or per-chromosome files), and
            {output_base}.haplotypes.index.duckdb.
        version_str: Version identifier embedded in sequence IDs (e.g. "1.0").
        merge: If True, include gnomAD single-variant sites above frequency_cutoff.
        frequency_cutoff: Minimum estimated gnomAD allele frequency for haplotype inclusion.
        split_contigs: If True, write one FASTA file per chromosome.
        tmp_dir: Temporary directory for Hail checkpoint files.
    """
    assert_directory_exists(haplotypes_table_path)
    assert_directory_exists(gnomad_va_file)
    assert_path_is_readable(reference_fasta)
    assert_fasta_indexed(reference_fasta)
    assert_path_is_writable(tmp_dir)
    assert_path_is_writable(output_base)

    hl.init(tmp_dir=tmp_dir)

    ht, pops_legend = build_haplotype_table(
        haplotypes_table_path=haplotypes_table_path,
        gnomad_va_file=gnomad_va_file,
        reference_fasta=reference_fasta,
        window_size=window_size,
        frequency_cutoff=frequency_cutoff,
        merge=merge,
        version_str=version_str,
        tmp_dir=tmp_dir,
    )

    file_suffix = ".haplotypes" if not merge else ".haplotypes_gnomad_merge"

    df = export_ht_to_dataframe(ht, output_base, file_suffix, pops_legend)
    write_fasta_files(df, output_base, file_suffix, split_contigs)
    create_duckdb_index(df, output_base, file_suffix, window_size, pops_legend, version_str)
