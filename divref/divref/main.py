import logging
import sys
from typing import Callable
from typing import List

import defopt

from divref.tools.extract_gnomad_afs import extract_gnomad_afs
from divref.tools.gnomad_hail_table_test_data import gnomad_hail_table_test_data

_tools: List[Callable[..., None]] = [
    extract_gnomad_afs,
    gnomad_hail_table_test_data,
]


def setup_logging(level: str = "INFO") -> None:
    """Set up basic logging to print to the console."""
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
