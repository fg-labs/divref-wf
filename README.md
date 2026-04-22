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
We checked all 'gnomAD_variant' variants on chr22 from DivRef 1.1 against:

- gnomAD 3.1.2 HGDP+1KG subset
- gnomAD 3.1.2 genomes (~76K genomes)
- gnomAD 4.1 joint (~730K exomes and ~76K genomes)

gnomAD 3.1.2 HGDP+1KG subset is the source used for the haplotypes present in DivRef 1.1.

```bash
pixi run -e analysis snakemake -j1 -s workflows/compare_divref_gnomad.smk
```

**DivRef 1.1 variants present in gnomAD datasets**

We found that all DivRef 1.1 'gnomAD_variant' variants were present in the gnomAD 3.1.2 HGDP+1KG subset and in the gnomAD 3.1.2 genomes set, while 16 were missing from the gnomAD 4.1 joint set.

We further compared the allele frequencies for the 5 populations recorded in the DivRef 1.1 DuckDB index against the frequencies for those populations in the two gnomAD sets, using the Hail tables as input.

- gnomAD 3.1.2 HGDP+1KG subset: within a rounding error of 5e06, for all variants.
- gnomAD 3.1.2 genomes: 1,400 variants with an AF difference >= 0.001, all of which had lower AF in the genomes set compared to HGDP+1KG, indicating more stringent filtering in the genomes set
- gnomAD 4.1 joint: 60,424 variants with an AF difference >= 0.001, evenly distributed both lower and higher AF

There are 506,983 DivRef 1.1 'gnomAD_variant' variants. When we run the original DivRef script [extract_gnomad_afs.py](https://github.com/e9genomics/human-diversity-reference/blob/main/scripts/extract_gnomad_afs.py) (uses gnomAD 3.1.2 HGDP+1KG sites as input, hard-coded) followed by lines 156-177 of [create_fasta_and_index.py](https://github.com/e9genomics/human-diversity-reference/blob/main/scripts/create_fasta_and_index.py) (lines 156-177) using the parameters specified in the [Makefile](https://github.com/e9genomics/human-diversity-reference/blob/main/Makefile), we also get 506,983 variants. 

We concluded that the DivRef 1.1 documentation was incorrect, and that the actual source of the gnomAD variants in the dataset was the gnomAD 3.1.2 HGDP+1KG subset, the same as for the haplotypes.
