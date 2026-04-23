####################################################################################################
# Generates a DivRef-format resource of human haplotypes.
#
# Final output is a set of per-chromosome FASTA files and a DuckDB index.
####################################################################################################

import os
from pathlib import Path
from snakemake.utils import validate

####################################################################################################
# Inputs
####################################################################################################


configfile: os.path.join(workflow.basedir, "config", "config.yml")


validate(config, os.path.join(workflow.basedir, "config", "config_schema.yml"))

VERSION: str = config["version"]

WORK_DIR: Path = Path(config["work_dir"])
TMP_DIR: Path = Path(config["tmp_dir"])

CHROMS: list[str] = config["chromosomes"]
POPS: list[str] = config["populations"]

REFERENCE_GENOME: str = config["reference_genome_base_name"]
REFERENCE_GENOME_URI: str = config["reference_genome_uri"]

# The HGDP+1KG phased BCF files are at
# "{HGDP_1KG_PHASED_BCF_PREFIX}.{chrom}.{HGDP_1KG_PHASED_BCF_SUFFIX}"
HGDP_1KG_PHASED_BCF_PREFIX: str = config["hgdp_1kg_phased_bcf_prefix"]
HGDP_1KG_PHASED_BCF_SUFFIX: str = config["hgdp_1kg_phased_bcf_suffix"]
HGDP_1KG_VARIANT_ANNOTATION_HAIL_TABLE: str = config["hgdp_1kg_variant_annotation_hail_table"]
HGDP_1KG_SAMPLE_METADATA_HAIL_TABLE: str = config["hgdp_1kg_sample_metadata_hail_table"]
HGDP_1KG_MIN_POP_AF_EXTRACT_GNOMAD_AFS: float = config["hgdp_1kg_min_pop_af_extract_gnomad_afs"]
HGDP_1KG_HAPLOTYPE_WINDOW_SIZE: int = config["hgdp_1kg_haplotype_window_size"]
HGDP_1KG_MIN_POP_AF_COMPUTE_HAPLOTYPES: float = config["hgdp_1kg_min_pop_af_compute_haplotypes"]
HGDP_1KG_MIN_EST_GNOMAD_HAPLOTYPE_AF: float = config["hgdp_1kg_min_estimated_gnomad_haplotype_af"]

SEQUENCE_WINDOW_SIZE: int = config["sequence_window_size"]

if HGDP_1KG_MIN_POP_AF_EXTRACT_GNOMAD_AFS > HGDP_1KG_MIN_POP_AF_COMPUTE_HAPLOTYPES:
    raise ValueError(
        f"hgdp_1kg_min_pop_af_extract_gnomad_afs ({HGDP_1KG_MIN_POP_AF_EXTRACT_GNOMAD_AFS}) "
        f"must be <= hgdp_1kg_min_pop_af_compute_haplotypes "
        f"({HGDP_1KG_MIN_POP_AF_COMPUTE_HAPLOTYPES})"
    )

VCF_EXTS: list[str] = [".vcf.gz", ".vcf.gz.tbi"]

####################################################################################################
# Rules
####################################################################################################


rule all:
    input:
        f"{WORK_DIR}/output/hgdp_1kg.haplotypes_gnomad_merge.index.duckdb",
        expand(
            f"{WORK_DIR}/output/hgdp_1kg.haplotypes_gnomad_merge.{{chrom}}.fasta",
            chrom=CHROMS,
        ),


####################################################################################################
# Extracts the phased genotypes for all HGDP+1KG samples in the specified locus.
#
# Removes the INFO field, since this is not required for haplotype computation (allele frequencies
# are re-annotated from the sites table), and it inflates the size on disk and subsequently the time
# for Hail to load and parse the VCF with the `divref compute-haplotypes` tool.
####################################################################################################
rule subset_phased_genotypes:
    output:
        vcf=f"{WORK_DIR}/inputs/hgdp_1kg.phased_genotypes.{{chrom}}.vcf.gz",
        tbi=f"{WORK_DIR}/inputs/hgdp_1kg.phased_genotypes.{{chrom}}.vcf.gz.tbi",
    log:
        "logs/generate_divref/subset_phased_genotypes.{chrom}.log",
    params:
        bcf=f"{HGDP_1KG_PHASED_BCF_PREFIX}{{chrom}}{HGDP_1KG_PHASED_BCF_SUFFIX}",
    shell:
        """
        (
            bcftools annotate \
                --remove INFO \
                --output-type z \
                --output {output.vcf} \
                --write-index=tbi \
                {params.bcf}
        ) &> {log}
        """


####################################################################################################
# Extracts allele frequencies from the HGDP+1KG gnomAD subset for the given populations and subsets
# to sites over the specified minimum allele frequency in at least one population.
####################################################################################################
rule extract_gnomad_afs:
    output:
        variant_ht=directory(f"{WORK_DIR}/inputs/hgdp_1kg.sites.{{chrom}}.ht"),
    log:
        "logs/generate_divref/extract_gnomad_afs.{chrom}.log",
    params:
        variant_ht=HGDP_1KG_VARIANT_ANNOTATION_HAIL_TABLE,
        freq_threshold=HGDP_1KG_MIN_POP_AF_EXTRACT_GNOMAD_AFS,
        populations=" ".join(POPS),
    shell:
        """
        (
            divref extract-gnomad-afs \
                --in-gnomad-sites-table {params.variant_ht} \
                --out-variant-annotation-table {output.variant_ht} \
                --contig {wildcards.chrom} \
                --freq-threshold {params.freq_threshold} \
                --populations {params.populations}
        ) &> {log}
        """


####################################################################################################
# Extracts selected fields from HGDP+1KG sample metadata.
####################################################################################################
rule extract_sample_metadata:
    output:
        sample_ht=directory(f"{WORK_DIR}/inputs/hgdp_1kg.sample_metadata.ht"),
    log:
        "logs/generate_divref/extract_sample_metadata.log",
    params:
        sample_ht=HGDP_1KG_SAMPLE_METADATA_HAIL_TABLE,
    shell:
        """
        (
            divref extract-sample-metadata \
                --in-gnomad-hgdp-sample-data {params.sample_ht} \
                --out-sample-metadata {output.sample_ht}
        ) &> {log}
        """


####################################################################################################
# Compute haplotypes from the HGDP+1KG filtered sites, sample metadata, and phased genotypes.
####################################################################################################
rule compute_haplotypes:
    input:
        vcf=f"{WORK_DIR}/inputs/hgdp_1kg.phased_genotypes.{{chrom}}.vcf.gz",
        tbi=f"{WORK_DIR}/inputs/hgdp_1kg.phased_genotypes.{{chrom}}.vcf.gz.tbi",
        variant_ht=f"{WORK_DIR}/inputs/hgdp_1kg.sites.{{chrom}}.ht",
        sample_ht=f"{WORK_DIR}/inputs/hgdp_1kg.sample_metadata.ht",
    output:
        haplotypes_ht=directory(f"{WORK_DIR}/haplotypes/hgdp_1kg.haplotypes.{{chrom}}.ht"),
    log:
        "logs/generate_divref/compute_haplotypes.{chrom}.log",
    params:
        window_size=HGDP_1KG_HAPLOTYPE_WINDOW_SIZE,
        variant_freq_threshold=HGDP_1KG_MIN_POP_AF_COMPUTE_HAPLOTYPES,
        haplotype_freq_threshold=HGDP_1KG_MIN_EST_GNOMAD_HAPLOTYPE_AF,
        output_base=f"{WORK_DIR}/haplotypes/hgdp_1kg.haplotypes.{{chrom}}",
    shell:
        """
        (
            divref compute-haplotypes \
                --vcfs-path {input.vcf} \
                --gnomad-va-file {input.variant_ht} \
                --gnomad-sa-file {input.sample_ht} \
                --window-size {params.window_size} \
                --variant-freq-threshold {params.variant_freq_threshold} \
                --haplotype-freq-threshold {params.haplotype_freq_threshold} \
                --output-base {params.output_base}
            
            # remove intermediate files
            rm -r {params.output_base}.[12].ht {params.output_base}.variants.ht
        ) &> {log}
        """


####################################################################################################
# Downloads and unzips the reference genome.
####################################################################################################
rule download_reference_genome:
    output:
        fasta=f"{WORK_DIR}/inputs/{REFERENCE_GENOME}.fasta",
    log:
        "logs/generate_divref/download_reference_genome.log",
    params:
        fasta_uri=REFERENCE_GENOME_URI,
    shell:
        """
        (
            gsutil -m cp {params.fasta_uri} {output.fasta}.gz
            gunzip {output.fasta}.gz
        ) &> {log}
        """


####################################################################################################
# Indexes the reference genome.
####################################################################################################
rule index_reference_genome:
    input:
        fasta=f"{WORK_DIR}/inputs/{REFERENCE_GENOME}.fasta",
    output:
        fai=f"{WORK_DIR}/inputs/{REFERENCE_GENOME}.fai",
    log:
        "logs/generate_divref/index_reference_genome.log",
    shell:
        """
        (
            samtools faidx \
                {input.fasta} \
                --output {output.fai}
        ) &> {log}
        """


####################################################################################################
# Writes a TSV listing per-chromosome haplotype and gnomAD sites Hail tables for the index builder.
####################################################################################################
rule create_table_pairs_tsv:
    input:
        haplotypes_hts=expand(
            f"{WORK_DIR}/haplotypes/hgdp_1kg.haplotypes.{{chrom}}.ht",
            chrom=CHROMS,
        ),
        sites_hts=expand(
            f"{WORK_DIR}/inputs/hgdp_1kg.sites.{{chrom}}.ht",
            chrom=CHROMS,
        ),
    output:
        tsv=f"{WORK_DIR}/inputs/hgdp_1kg.table_pairs.tsv",
    run:
        with open(output.tsv, "w") as f:
            f.write("contig\thaplotype_table_path\tsites_table_path\n")
            for chrom in CHROMS:
                haplotype_ht = f"{WORK_DIR}/haplotypes/hgdp_1kg.haplotypes.{chrom}.ht"
                sites_ht = f"{WORK_DIR}/inputs/hgdp_1kg.sites.{chrom}.ht"
                f.write(f"{chrom}\t{haplotype_ht}\t{sites_ht}\n")


####################################################################################################
# Build the DivRef DuckDB index from all per-chromosome haplotype and gnomAD sites Hail tables.
####################################################################################################
rule create_divref_index:
    input:
        table_pairs_tsv=f"{WORK_DIR}/inputs/hgdp_1kg.table_pairs.tsv",
        fasta=f"{WORK_DIR}/inputs/{REFERENCE_GENOME}.fasta",
        fai=f"{WORK_DIR}/inputs/{REFERENCE_GENOME}.fai",
    output:
        duckdb=f"{WORK_DIR}/output/hgdp_1kg.haplotypes_gnomad_merge.index.duckdb",
    log:
        "logs/generate_divref/create_divref_index.log",
    params:
        window_size=SEQUENCE_WINDOW_SIZE,
        output_base=f"{WORK_DIR}/output/hgdp_1kg",
        version=VERSION,
        tmp_dir=TMP_DIR,
    shell:
        """
        (
            divref create-duckdb-index \
                --in-table-pairs-tsv {input.table_pairs_tsv} \
                --reference-fasta {input.fasta} \
                --window-size {params.window_size} \
                --output-base {params.output_base} \
                --version {params.version} \
                --tmp-dir {params.tmp_dir}
        ) &> {log}
        """


####################################################################################################
# Write per-chromosome FASTA files from the DivRef DuckDB index.
####################################################################################################
rule create_divref_fasta:
    input:
        duckdb=f"{WORK_DIR}/output/hgdp_1kg.haplotypes_gnomad_merge.index.duckdb",
    output:
        fastas=expand(
            f"{WORK_DIR}/output/hgdp_1kg.haplotypes_gnomad_merge.{{chrom}}.fasta",
            chrom=CHROMS,
        ),
    log:
        "logs/generate_divref/create_divref_fasta.log",
    params:
        output_base=f"{WORK_DIR}/output/hgdp_1kg.haplotypes_gnomad_merge",
    shell:
        """
        (
            divref create-divref-fasta \
                --duckdb-path {input.duckdb} \
                --output-base {params.output_base}
        ) &> {log}
        """
