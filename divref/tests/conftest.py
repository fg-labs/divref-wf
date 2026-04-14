"""Shared pytest fixtures for the divref test suite."""

from collections.abc import Generator
from pathlib import Path

import hail as hl
import pytest


@pytest.fixture(scope="session")
def hail_context() -> Generator[None, None, None]:
    """
    Initialize a local Hail context for the test session.

    Starts a local Spark/JVM context once per session and stops it on exit.
    Tests that use Hail expressions must request this fixture explicitly.
    Requires Java 11+ to be available in PATH (provided by pixi or setup-java).
    """
    hl.init(quiet=True)
    yield
    hl.stop()


@pytest.fixture
def datadir() -> Path:
    """Path to tests/data."""
    return Path(__file__).parent / "data"


@pytest.fixture
def hail_reference_genome(hail_context: None) -> hl.ReferenceGenome:  # noqa: ARG001
    """A small custom reference genome for use in testing."""
    contigs: list[str] = ["chr1"]
    lengths: dict[str, int] = {"chr1": 1000}

    reference_genome = hl.ReferenceGenome(
        "test_chr1", contigs, lengths, x_contigs=[], y_contigs=[], mt_contigs=[]
    )
    return reference_genome
