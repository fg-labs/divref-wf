"""Tool to write per-chromosome FASTA files from a DivRef DuckDB index."""

import logging
from pathlib import Path

import duckdb
import polars
from fgpyo.io import assert_path_is_readable
from fgpyo.io import assert_path_is_writable

logger = logging.getLogger(__name__)


def _write_fasta_files(df: polars.DataFrame, output_base: Path) -> None:
    """
    Write one FASTA file per chromosome to {output_base}.{chrom}.fasta.

    Args:
        df: DataFrame with sequence_id, sequence, and variants columns.
        output_base: Base path; chromosome name is appended as a suffix.
    """
    for chrom in sorted(df["contig"].unique().to_list()):
        logger.info("Creating FASTA for chromosome %s", chrom)
        df_chrom = df.filter(df["contig"] == chrom)
        out_path = Path(f"{output_base}.{chrom}.fasta")
        with open(out_path, "w") as fasta_out:
            for sequence_id, sequence in df_chrom.select("sequence_id", "sequence").iter_rows():
                fasta_out.write(f">{sequence_id}\n{sequence}\n")


def create_divref_fasta(
    *,
    duckdb_path: Path,
    output_base: Path,
) -> None:
    """
    Write per-chromosome FASTA files from a DivRef DuckDB index.

    Reads sequence_id, sequence, and variants from the sequences table and writes one FASTA
    file per chromosome to {output_base}.{chrom}.fasta.

    Args:
        duckdb_path: Path to an existing DivRef DuckDB index.
        output_base: Base path for output FASTA files; chromosome name is appended as a suffix.
    """
    assert_path_is_readable(duckdb_path)
    assert_path_is_writable(output_base)
    con = duckdb.connect(str(duckdb_path), read_only=True)
    df = con.execute("SELECT sequence_id, sequence, contig FROM sequences").pl()
    con.close()
    _write_fasta_files(df, output_base)
