"""Tool to build DivRef FASTA sequences and DuckDB index from haplotype Hail tables."""

import logging
import os

import duckdb
import hail as hl
import polars

from divref.alias import HailPath
from divref.haplotype import get_haplo_sequence
from divref.haplotype import split_haplotypes

logger = logging.getLogger(__name__)


def create_fasta_and_index(
    *,
    haplotypes_table_path: HailPath,
    gnomad_va_file: HailPath,
    reference_fasta: HailPath,
    window_size: int,
    output_base: HailPath,
    version_str: str,
    merge: bool = False,
    frequency_cutoff: float = 0.005,
    split_contigs: bool = False,
    tmp_dir: HailPath = "/tmp",
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
    hl.init(tmp_dir=tmp_dir)

    ht = hl.read_table(haplotypes_table_path).key_by()
    va = hl.read_table(gnomad_va_file)
    pops_legend: list[str] = va.pops.collect()[0]

    hl.get_reference("GRCh38").add_sequence(reference_fasta)

    logger.info(
        "Haplotype table contains %d unique haplotypes above frequency threshold", ht.count()
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
    ).export(output_base + f"{file_suffix}.tsv.bgz")

    df = polars.read_csv(
        output_base + f"{file_suffix}.tsv.bgz",
        separator="\t",
        schema_overrides={"sequence_id": polars.String},
    )

    if split_contigs:
        df = df.with_columns(contig=df["variants"].str.split(":").list.get(0))
        for chrom in df["contig"].unique().to_list():
            logger.info("Creating FASTA for chromosome %s", chrom)
            df2 = df.filter(df["contig"] == chrom)
            with open(output_base + f"{file_suffix}.{chrom}.fasta", "w") as fasta_out:
                for sequence, sequence_id in df2.select("sequence", "sequence_id").iter_rows():
                    fasta_out.write(f">{sequence_id}\n{sequence}\n")
    else:
        logger.info("Creating FASTA")
        with open(output_base + f"{file_suffix}.fasta", "w") as fasta_out:
            for sequence, sequence_id in df.select("sequence", "sequence_id").iter_rows():
                fasta_out.write(f">{sequence_id}\n{sequence}\n")

    duckdb_file = output_base + f"{file_suffix}.index.duckdb"
    if os.path.exists(duckdb_file):
        os.remove(duckdb_file)
    con = duckdb.connect(duckdb_file)
    con.execute("CREATE TABLE sequences AS SELECT * FROM df")
    con.execute("CREATE INDEX idx_sequence_id ON sequences(sequence_id)")
    con.execute(f"CREATE TABLE window_size AS SELECT {window_size} AS window_size")
    con.execute(f"CREATE TABLE pops_legend AS SELECT {pops_legend} AS pops_legend")
    con.execute(f"CREATE TABLE VERSION AS SELECT '{version_str}' AS version")
    con.close()
