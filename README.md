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

The environment and dependencies are automatically created and installed when calling `pixi install` for the first time.

To enable access to Hail tables via the GCS Connector, run `pixi run setup-gcs`.

You will need to log in to GCS before running any of the Hail-dependent tools.

```bash
gcloud auth application-default login
```

## Source Data

### gnomAD 3.1.2 HGDP+1KG individual-level genotypes and sample metadata

[Data description](https://gnomad.broadinstitute.org/news/2021-10-gnomad-v3-1-2-minor-release/)

## Analysis

### Additional Environment Requirements

To install R packages not available as conda-forge builds for all platforms (duckdb, duckplyr), run

```bash
pixi run -e analysis setup-r-packages
```

### Compare DivRef 1.1 against different gnomAD releases

[DivRef 1.1](https://zenodo.org/records/14802613) states that:

> DivRef is constructed by computing empirical phased haplotypes within 25 BPs over 0.5% allele frequency from the Human Genome Diversity Panel (HGDP) using the phased Hail dataset provided by the gnomAD team at the Broad Institute, merged with single variants over 0.5% AF from the gnomAD v4.1.0 summary release.

Some gnomAD v4.1.0 variants that we expected to see represented in DivRef 1.1 were not.
We checked all 'gnomAD_variant' variants on chr22 from DivRef 1.1 against the gnomAD 4.1 joint (exomes and genomes) variants and the gnomAD 3.1.2 HGDP+1KG subset variants.
The latter is the source used for the haplotypes present in DivRef.

```bash
pixi run -e analysis snakemake -j1 -s workflows/compare_divref_gnomad.smk
```

We found that all DivRef 1.1 variants were present in the gnomAD 3.1.2 HGDP+1KG subset, while 16 were missing from the gnomAD 4.1 joint set.
We further compared the allele frequencies for the 5 populations recorded in the DivRef 1.1 DuckDB index against the frequencies for those populations in the two gnomAD sets, using the Hail tables as input.

The allele frequencies between DivRef 1.1 and the gnomAD 3.1.2 HGDP+1KG subset were very close, within a rounding error of 5e06, for all variants.

The allele frequencies between DivRef 1.1 and gnomAD 4.1 joint set were significantly different. 60,424 variants were found with an allele frequence difference >= 0.001.

We concluded that the DivRef 1.1 documentation was incorrect, and that the actual source of the gnomAD variants in the dataset was the gnomAD 3.1.2 HGDP+1KG subset, the same as for the haplotypes.
