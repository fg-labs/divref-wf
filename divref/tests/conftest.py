"""Shared pytest fixtures for the divref test suite."""

from collections.abc import Generator

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
