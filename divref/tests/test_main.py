from collections.abc import Callable
from inspect import isfunction

import pytest
from defopt import signature

from divref import main


@pytest.mark.parametrize("tool", main._tools)
def test_tools_are_defined(tool: Callable[..., None]) -> None:
    """Test that all command line tools passed to defopt are defined functions."""
    assert isfunction(tool)


@pytest.mark.parametrize("tool", main._tools)
def test_tools_have_valid_docstrings(tool: Callable[..., None]) -> None:
    """Test that all command line tools have a valid defopt docstring."""
    try:
        signature(tool)
    except TypeError:
        raise AssertionError(f"defopt could not parse docstring for {tool.__name__}") from None


@pytest.mark.parametrize("tool", main._tools, ids=lambda t: t.__name__)
def test_cli_help_builds_every_subcommand_parser(
    tool: Callable[..., None], monkeypatch: pytest.MonkeyPatch
) -> None:
    """
    `--help` builds each subcommand's argparse parser.

    defopt builds every subcommand parser eagerly, so this fails at parser-build time if any
    parameter type is unresolved.
    """
    subcommand = tool.__name__.replace("_", "-")
    monkeypatch.setattr("sys.argv", ["divref", subcommand, "--help"])
    with pytest.raises(SystemExit) as excinfo:
        main.run()
    assert excinfo.value.code == 0
