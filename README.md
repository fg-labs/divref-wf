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

1. Download the DivRef 1.1 index file.

```bash
mkdir -p data/analysis/input
wget -O data/analysis/input/DivRef-v1.1.haplotypes_gnomad_merge.index.duckdb \
  https://zenodo.org/records/14802613/files/DivRef-v1.1.haplotypes_gnomad_merge.index.duckdb
```

2. Extract locus, alleles, and default DivRef population AFs from the gnomAD 4.1 joint (exomes and genomes) call set and the gnomAD 3.1.2 HGDP+1KG subset call set. Restrict to `chr22`.

```bash
pixi run divref extract-gnomad-single-afs \
  --contig chr22 \
  --gnomad-version JOINT_41 \
  --out-sites-hail-table data/analysis/input/chr22.joint_41.ht \
  --out-sites-tsv data/analysis/input/chr22.joint_41.tsv

pixi run divref extract-gnomad-single-afs \
  --contig chr22 \
  --gnomad-version HGDP_1KG_312 \
  --out-sites-hail-table data/analysis/input/chr22.hgdp_1kg_312.ht \
  --out-sites-tsv data/analysis/input/chr22.hgdp_1kg_312.ts
```

3. Compare DivRef 1.1 chr22 gnomAD single variant sites to each of the two gnomAD releases.

```bash
mkdir -p data/analysis/compare_divref_gnomad

pixi run Rscript scripts/compare_divref_gnomad.R \
  --contig chr22 \
  --divref_duckdb data/analysis/input/DivRef-v1.1.haplotypes_gnomad_merge.index.duckdb \
  --gnomad_tsv data/analysis/input/chr22.joint_41.tsv \
  --gnomad_label "gnomAD 4.1 joint AF" \
  --output_base data/analysis/compare_divref_gnomad/chr22.joint_41

pixi run Rscript scripts/compare_divref_gnomad.R \
  --contig chr22 \
  --divref_duckdb data/analysis/input/DivRef-v1.1.haplotypes_gnomad_merge.index.duckdb \
  --gnomad_tsv data/analysis/input/chr22.hgdp_1kg_312.tsv \
  --gnomad_label "gnomAD 3.1.2 HGDP+1KG AF" \
  --output_base data/analysis/compare_divref_gnomad/chr22.hgdp_1kg_312
```
