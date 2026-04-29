"""Tool to write per-chromosome FASTA files from a DivRef DuckDB index."""

import logging
from pathlib import Path

import duckdb
from fgpyo.io import assert_path_is_readable
from fgpyo.io import assert_path_is_writable

logger = logging.getLogger(__name__)


def create_divref_fasta(*, duckdb_path: Path, output_base: Path, contigs: list[str]) -> None:
    """
    Write per-chromosome FASTA files from a DivRef DuckDB index.

    Reads sequence_id, sequence, and variants from the sequences table and writes one FASTA
    file per chromosome to {output_base}.{chrom}.fasta.

    Args:
        duckdb_path: Path to an existing DivRef DuckDB index.
        output_base: Base path for output FASTA files; chromosome name is appended as a suffix.
        contigs: List of contigs to write FASTA files for.
    """
    assert_path_is_readable(duckdb_path)
    assert_path_is_writable(output_base)

    if len(contigs) == 0:
        raise ValueError("Contig list must be provided.")

    with duckdb.connect(str(duckdb_path), read_only=True) as con:
        for contig in contigs:
            df = con.execute(
                "SELECT sequence_id, sequence FROM sequences WHERE contig = $contig",
                {"contig": contig},
            ).pl()

            out_path = Path(f"{output_base}.{contig}.fasta")
            assert_path_is_writable(out_path)
            logger.info(f"Creating FASTA for chromosome {contig} at {out_path}")
            with out_path.open("w") as fh:
                for sequence_id, sequence in df.select("sequence_id", "sequence").iter_rows():
                    fh.write(f">{sequence_id}\n{sequence}\n")
