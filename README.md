[![CI](https://github.com/fg-labs/divref-wf/actions/workflows/python_package.yml/badge.svg?branch=main)](https://github.com/fg-labs/divref-wf/actions/workflows/python_package.yml?query=branch%3Amain)
[![Python Versions](https://img.shields.io/badge/python-3.12_|_3.13-blue)](https://github.com/fg-labs/divref-wf)
[![MyPy Checked](http://www.mypy-lang.org/static/mypy_badge.svg)](http://mypy-lang.org/)
[![uv](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/uv/main/assets/badge/v0.json)](https://github.com/astral-sh/uv)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://docs.astral.sh/ruff/)
[![Pixi](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/prefix-dev/pixi/main/assets/badge/v0.json)](https://pixi.sh)

# Snakemake workflow implementation to create DivRef-style resource

This workflow is inspired by the [DivRef](https://github.com/e9genomics/human-diversity-reference) repository which is used to generate a bundle of FASTA sequences and a corresponding DuckDB index of common human variation.

The original implementation is via a set of standalone Python scripts and a Makefile.
This implementation:
1. Wraps the Python scripts in a toolkit with added typing, improved parameterization, and added unit testing.
2. Adds a Snakemake workflow and associated configuration to drive the resource generation process.

## Set up Environment

The environment for this analysis is managed using `pixi`.
Follow the developer [instructions](https://pixi.sh/latest/installation/) to install `pixi`.

The environment and dependencies are automatically created and installed when calling `pixi run <cmd>` for the first time.
