import logging
import sys
from typing import Callable
from typing import List

import defopt

from divref.tools.compute_haplotype_statistics import compute_haplotype_statistics
from divref.tools.compute_haplotypes import compute_haplotypes
from divref.tools.compute_variation_ratios import compute_variation_ratios
from divref.tools.create_fasta_and_index import create_fasta_and_index
from divref.tools.create_gnomad_sites_vcf import create_gnomad_sites_vcf
from divref.tools.extract_gnomad_afs import extract_gnomad_afs
from divref.tools.extract_sample_metadata import extract_sample_metadata
from divref.tools.gnomad_hail_table_test_data import gnomad_hail_table_test_data
from divref.tools.remap_divref import remap_divref

_tools: List[Callable[..., None]] = [
    compute_haplotypes,
    compute_haplotype_statistics,
    compute_variation_ratios,
    create_fasta_and_index,
    create_gnomad_sites_vcf,
    extract_gnomad_afs,
    extract_sample_metadata,
    gnomad_hail_table_test_data,
    remap_divref,
]


def setup_logging(level: str = "INFO") -> None:
    """
    Set up basic logging to print to the console.

    Args:
        level: Logging level string (e.g. "INFO", "DEBUG", "WARNING").
    """
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(name)s:%(funcName)s:%(lineno)s [%(levelname)s]: %(message)s",
    )


def run() -> None:
    """Set up logging, then hand over to defopt for running command line tools."""
    setup_logging()
    logger = logging.getLogger("divref")
    logger.info("Executing: " + " ".join(sys.argv))
    defopt.run(
        funcs=_tools,
        argv=sys.argv[1:],
    )
    logger.info("Finished executing successfully.")
