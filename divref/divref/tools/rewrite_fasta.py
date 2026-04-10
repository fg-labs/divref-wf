"""Tool to filter a FASTA file to canonical chromosomes only."""

from pathlib import Path

import tqdm


def rewrite_fasta(*, fasta_path: Path, output_path: Path) -> None:
    """
    Rewrite a FASTA file keeping only canonical chromosomes (chr1-22, X, Y, MT).

    Filters out alt contigs, decoy sequences, and any other non-canonical sequences.
    Reads the input file line by line so it can handle arbitrarily large FASTA files.

    Args:
        fasta_path: Path to the input FASTA file.
        output_path: Path to write the filtered FASTA file.
    """
    contigs_to_keep = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrMT"}

    keep = False
    with open(fasta_path) as f, open(output_path, "w") as out:
        for line in tqdm.tqdm(f):
            if line[0] == ">":
                contig = line.split()[0][1:]
                keep = contig in contigs_to_keep
                if keep:
                    out.write(line)
            elif keep:
                out.write(line)
