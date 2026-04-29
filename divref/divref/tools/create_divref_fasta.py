"""Tool to write per-chromosome FASTA files from a DivRef DuckDB index."""

import logging
from collections.abc import Iterator
from pathlib import Path

import duckdb
import polars
from fgpyo.io import assert_path_is_readable
from fgpyo.io import assert_path_is_writable

logger = logging.getLogger(__name__)


def create_divref_fasta(
    *,
    duckdb_path: Path,
    output_base: Path,
    contigs: list[str],
    polars_chunk_size: int = 100_000,
) -> None:
    """
    Write per-chromosome FASTA files from a DivRef DuckDB index.

    Streams sequences for each contig out of the DuckDB sequences table in batches of
    ``polars_chunk_size`` rows and writes them directly to ``{output_base}.{contig}.fasta`` without
    materialising the per-contig result as a single in-process DataFrame. If a requested contig
    has no rows in the index, an empty FASTA is written and a warning is logged.

    Args:
        duckdb_path: Path to an existing DivRef DuckDB index.
        output_base: Base path for output FASTA files; chromosome name is appended as a suffix.
        contigs: List of contigs to write FASTA files for.
        polars_chunk_size: Maximum number of rows per polars DataFrame batch read from DuckDB.
    """
    assert_path_is_readable(duckdb_path)
    assert_path_is_writable(output_base)

    if len(contigs) == 0:
        raise ValueError("Contig list must be provided.")

    with duckdb.connect(str(duckdb_path), read_only=True) as conn:
        for contig in contigs:
            out_path = Path(f"{output_base}.{contig}.fasta")
            assert_path_is_writable(out_path)
            logger.info(f"Creating FASTA for chromosome {contig} at {out_path}")
            rows_written: int = 0
            with out_path.open("w") as fh:
                for df in iter_sequence_chunks(
                    conn=conn, contig=contig, chunk_size=polars_chunk_size
                ):
                    for sequence_id, sequence in df.iter_rows():
                        fh.write(f">{sequence_id}\n{sequence}\n")
                    rows_written += df.height
            if rows_written == 0:
                logger.warning(
                    f"No sequences found for contig {contig}; wrote empty FASTA at {out_path}"
                )
            else:
                logger.info(f"Wrote {rows_written} sequences for {contig}")


def iter_sequence_chunks(
    *,
    conn: duckdb.DuckDBPyConnection,
    contig: str,
    chunk_size: int,
) -> Iterator[polars.DataFrame]:
    """
    Yield polars DataFrames of ``(sequence_id, sequence)`` rows for one contig.

    Streams the DuckDB result set as Arrow record batches and converts each batch to polars,
    bounding in-process memory by ``chunk_size`` rows.

    Args:
        conn: Open DuckDB connection to the DivRef index.
        contig: Contig to filter on.
        chunk_size: Maximum rows per yielded DataFrame.

    Yields:
        Polars DataFrame batches with ``sequence_id`` and ``sequence`` columns.
    """
    result = conn.execute(
        "SELECT sequence_id, sequence FROM sequences WHERE contig = $contig",
        {"contig": contig},
    )
    for batch in result.fetch_record_batch(chunk_size):
        df = polars.from_arrow(batch)
        if isinstance(df, polars.DataFrame) and df.height > 0:
            yield df
