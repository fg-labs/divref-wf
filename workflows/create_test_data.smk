####################################################################################################
# Creates the test data in divref/tests/data.
####################################################################################################

from pathlib import Path

####################################################################################################
# Inputs
####################################################################################################

OUTPUT_DIR: Path = Path("divref/tests/data")
LOCUS_CHROM: str = "chr1"
LOCUS: str = "chr1:100001-200000"
LOCUS_FILENAME: str = "chr1_100001_200000"
# GRCh38 non-PAR boundaries for the two sex-chromosome loci below (each window sits well inside
# non-PAR): PAR1 ends at 2,781,479; chrX PAR2 starts at 155,701,383; chrY PAR2 at 56,887,903.
# chrX non-PAR locus exercises the haploid-male ploidy correction in `compute_haplotypes`.
CHRX_LOCUS_CHROM: str = "chrX"
CHRX_LOCUS: str = "chrX:50000000-50025000"
CHRX_LOCUS_FILENAME: str = "chrX_50000000_50025000"
CHRX_BCF_NONPAR: str = (
    "gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/"
    "hgdp1kgp_chrX_non_par.full.shapeit5_rare.bcf"
)
# chrY non-PAR locus exercises haplotype construction from unphased male-only genotypes.
# 27 PASS+common sites fall in-window; 94 of the 480 fixture samples (all male, since only males
# have called chrY genotypes) carry >=2 alt variants there, so haplotypes form.
CHRY_LOCUS_CHROM: str = "chrY"
CHRY_LOCUS: str = "chrY:2900000-2925000"
CHRY_LOCUS_FILENAME: str = "chrY_2900000_2925000"
CHRY_VCF: str = (
    "gs://gcp-public-data--gnomad/release/3.1.2/vcf/genomes/"
    "gnomad.genomes.v3.1.2.hgdp_tgp.chrY.vcf.bgz"
)
MIN_POP_AF_EXTRACT_GNOMAD_AFS: float = 0.001
MIN_POP_AF_COMPUTE_HAPLOTYPES: float = 0.005
WINDOW_SIZE_COMPUTE_HAPLOTYPES: int = 25
# Hail-using divref tools require a GCS credentials path when reading from local-only Hail
# tables, because hail_init currently sets `use_s3=False` by default and asserts the path is
# present. Passed explicitly to the rules below (the first subset rule relies on the tool's
# matching default of this same path).
GCS_CREDENTIALS_PATH: str = "~/.config/gcloud/application_default_credentials.json"

####################################################################################################
# Rules
####################################################################################################


rule all:
    input:
        f"{OUTPUT_DIR}/{LOCUS_FILENAME}.ht",
        f"{OUTPUT_DIR}/hgdp_1kg_sample_metadata.ht",
        f"{OUTPUT_DIR}/samples.txt",
        f"{OUTPUT_DIR}/{LOCUS_FILENAME}.vcf.gz",
        f"{OUTPUT_DIR}/{LOCUS_FILENAME}.vcf.gz.tbi",
        f"{OUTPUT_DIR}/{LOCUS_FILENAME}.gnomad_afs.ht",
        f"{OUTPUT_DIR}/hgdp_1kg_sample_metadata.extract.ht",
        f"{OUTPUT_DIR}/{LOCUS_FILENAME}_haplotypes.ht",
        f"{OUTPUT_DIR}/{CHRX_LOCUS_FILENAME}.ht",
        f"{OUTPUT_DIR}/{CHRX_LOCUS_FILENAME}.vcf.gz",
        f"{OUTPUT_DIR}/{CHRX_LOCUS_FILENAME}.vcf.gz.tbi",
        f"{OUTPUT_DIR}/{CHRX_LOCUS_FILENAME}.gnomad_afs.ht",
        f"{OUTPUT_DIR}/{CHRY_LOCUS_FILENAME}.ht",
        f"{OUTPUT_DIR}/{CHRY_LOCUS_FILENAME}.vcf.gz",
        f"{OUTPUT_DIR}/{CHRY_LOCUS_FILENAME}.vcf.gz.tbi",
        f"{OUTPUT_DIR}/{CHRY_LOCUS_FILENAME}.gnomad_afs.ht",


####################################################################################################
# Extracts all gnomAD HGDP+1KG variants in the specified locus from
# gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_variant_annotations.ht.
#
# Extracts selected fields from gnomAD sample metadata required by downstream tools from
# gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.ht.
####################################################################################################
rule subset_gnomad_hail_tables:
    output:
        chr1_variant_ht=directory(f"{OUTPUT_DIR}/{LOCUS_FILENAME}.ht"),
        chrx_variant_ht=directory(f"{OUTPUT_DIR}/{CHRX_LOCUS_FILENAME}.ht"),
        chry_variant_ht=directory(f"{OUTPUT_DIR}/{CHRY_LOCUS_FILENAME}.ht"),
        sample_ht=directory(f"{OUTPUT_DIR}/hgdp_1kg_sample_metadata.ht"),
        samples_txt=f"{OUTPUT_DIR}/samples.txt",
    log:
        "logs/create_test_data/subset_gnomad_hail_tables.log",
    params:
        chr1_locus=LOCUS,
        chrx_locus=CHRX_LOCUS,
        chry_locus=CHRY_LOCUS,
        subsample_seed=0,
    shell:
        """
        (
            divref gnomad-hail-table-test-data \
                --loci {params.chr1_locus} {params.chrx_locus} {params.chry_locus} \
                --out-variant-annotation-dir {OUTPUT_DIR} \
                --out-sample-metadata {output.sample_ht} \
                --out-samples-txt {output.samples_txt} \
                --subsample-seed {params.subsample_seed}
        ) &> {log}
        """


####################################################################################################
# Extracts the phased genotypes for all HGDP+1KG samples in the specified locus.
####################################################################################################
rule subset_phased_genotypes:
    input:
        samples_txt=f"{OUTPUT_DIR}/samples.txt",
    output:
        vcf=f"{OUTPUT_DIR}/{LOCUS_FILENAME}.vcf.gz",
        tbi=f"{OUTPUT_DIR}/{LOCUS_FILENAME}.vcf.gz.tbi",
    log:
        f"logs/create_test_data/subset_phased_genotypes.{LOCUS_FILENAME}.log",
    params:
        locus=LOCUS,
        bcf=f"gs://gcp-public-data--gnomad/resources/hgdp_1kg/phased_haplotypes_v2/hgdp1kgp_{LOCUS_CHROM}.filtered.SNV_INDEL.phased.shapeit5.bcf",
    shell:
        """
        (
            bcftools view \
                --regions {params.locus} \
                --samples-file {input.samples_txt} \
                --force-samples \
                --output-type z \
                --output {output.vcf} \
                --write-index=tbi \
                {params.bcf}
        ) &> {log}
        """


####################################################################################################
# Extracts allele frequencies for the default populations and subsets to sites over the specified
# minimum frequency.
####################################################################################################
rule extract_gnomad_afs:
    input:
        variant_ht=f"{OUTPUT_DIR}/{LOCUS_FILENAME}.ht",
    output:
        variant_ht=directory(f"{OUTPUT_DIR}/{LOCUS_FILENAME}.gnomad_afs.ht"),
    log:
        f"logs/create_test_data/extract_gnomad_afs.{LOCUS_FILENAME}.log",
    params:
        contig=LOCUS_CHROM,
        freq_threshold=MIN_POP_AF_EXTRACT_GNOMAD_AFS,
        gcs_credentials_path=GCS_CREDENTIALS_PATH,
    shell:
        """
        (
            divref extract-gnomad-afs \
                --in-gnomad-sites-table {input.variant_ht} \
                --out-variant-annotation-table {output.variant_ht} \
                --contig {params.contig} \
                --freq-threshold {params.freq_threshold} \
                --gcs-credentials-path {params.gcs_credentials_path}
        ) &> {log}
        """


####################################################################################################
# Extracts selected fields from sample metadata.
####################################################################################################
rule extract_sample_metadata:
    input:
        sample_ht=f"{OUTPUT_DIR}/hgdp_1kg_sample_metadata.ht",
    output:
        sample_ht=directory(f"{OUTPUT_DIR}/hgdp_1kg_sample_metadata.extract.ht"),
    log:
        f"logs/create_test_data/extract_sample_metadata.{LOCUS_FILENAME}.log",
    params:
        gcs_credentials_path=GCS_CREDENTIALS_PATH,
    shell:
        """
        (
            divref extract-sample-metadata \
                --in-gnomad-hgdp-sample-data {input.sample_ht} \
                --out-sample-metadata {output.sample_ht} \
                --gcs-credentials-path {params.gcs_credentials_path}
        ) &> {log}
        """


####################################################################################################
# Compute haplotypes from the sites and phased genotypes.
####################################################################################################
rule compute_haplotypes:
    input:
        vcf=f"{OUTPUT_DIR}/{LOCUS_FILENAME}.vcf.gz",
        tbi=f"{OUTPUT_DIR}/{LOCUS_FILENAME}.vcf.gz.tbi",
        variant_ht=f"{OUTPUT_DIR}/{LOCUS_FILENAME}.gnomad_afs.ht",
        sample_ht=f"{OUTPUT_DIR}/hgdp_1kg_sample_metadata.extract.ht",
    output:
        haplotypes_ht=directory(f"{OUTPUT_DIR}/{LOCUS_FILENAME}_haplotypes.ht"),
    log:
        f"logs/create_test_data/compute_haplotypes.{LOCUS_FILENAME}.log",
    params:
        window_size=WINDOW_SIZE_COMPUTE_HAPLOTYPES,
        freq_threshold=MIN_POP_AF_COMPUTE_HAPLOTYPES,
        output_base=f"{OUTPUT_DIR}/{LOCUS_FILENAME}_haplotypes",
    shell:
        """
        (
            divref compute-haplotypes \
                --vcfs-path {input.vcf} \
                --gnomad-va-file {input.variant_ht} \
                --gnomad-sa-file {input.sample_ht} \
                --window-size {params.window_size} \
                --variant-freq-threshold {params.freq_threshold} \
                --haplotype-freq-threshold {params.freq_threshold} \
                --output-base {params.output_base}
        ) &> {log}
        """


####################################################################################################
# Extracts phased genotypes for the chrX non-PAR test locus directly from the non-PAR BCF.
####################################################################################################
rule subset_phased_genotypes_chrX:
    input:
        samples_txt=f"{OUTPUT_DIR}/samples.txt",
    output:
        vcf=f"{OUTPUT_DIR}/{CHRX_LOCUS_FILENAME}.vcf.gz",
        tbi=f"{OUTPUT_DIR}/{CHRX_LOCUS_FILENAME}.vcf.gz.tbi",
    log:
        f"logs/create_test_data/subset_phased_genotypes.{CHRX_LOCUS_FILENAME}.log",
    params:
        locus=CHRX_LOCUS,
        bcf=CHRX_BCF_NONPAR,
    shell:
        """
        (
            bcftools view \
                --regions {params.locus} \
                --samples-file {input.samples_txt} \
                --force-samples \
                --output-type z \
                --output {output.vcf} \
                --write-index=tbi \
                {params.bcf}
        ) &> {log}
        """


####################################################################################################
# Extract gnomAD allele frequencies for the chrX non-PAR test locus.
####################################################################################################
rule extract_gnomad_afs_chrX:
    input:
        variant_ht=f"{OUTPUT_DIR}/{CHRX_LOCUS_FILENAME}.ht",
    output:
        variant_ht=directory(f"{OUTPUT_DIR}/{CHRX_LOCUS_FILENAME}.gnomad_afs.ht"),
    log:
        f"logs/create_test_data/extract_gnomad_afs.{CHRX_LOCUS_FILENAME}.log",
    params:
        contig=CHRX_LOCUS_CHROM,
        freq_threshold=MIN_POP_AF_EXTRACT_GNOMAD_AFS,
        gcs_credentials_path=GCS_CREDENTIALS_PATH,
    shell:
        """
        (
            divref extract-gnomad-afs \
                --in-gnomad-sites-table {input.variant_ht} \
                --out-variant-annotation-table {output.variant_ht} \
                --contig {params.contig} \
                --freq-threshold {params.freq_threshold} \
                --gcs-credentials-path {params.gcs_credentials_path}
        ) &> {log}
        """


####################################################################################################
# Extracts genotypes for the chrY non-PAR test locus from the raw gnomAD chrY VCF.
#
# Unlike the chr1/chrX BCFs (already PASS-only and phased), this is a raw multi-sample callset:
# it must be filtered to PASS sites and stripped of INFO/extra FORMAT fields before use.
# The chrY VCF carries all samples; females appear as `./.` missing calls, not absent, since only
# males have called chrY genotypes. `--force-samples` is defensive here (it tolerates any sample in
# `samples.txt` not found in the VCF), not required to drop females.
####################################################################################################
rule subset_genotypes_chrY:
    input:
        samples_txt=f"{OUTPUT_DIR}/samples.txt",
    output:
        vcf=f"{OUTPUT_DIR}/{CHRY_LOCUS_FILENAME}.vcf.gz",
        tbi=f"{OUTPUT_DIR}/{CHRY_LOCUS_FILENAME}.vcf.gz.tbi",
    log:
        f"logs/create_test_data/subset_genotypes.{CHRY_LOCUS_FILENAME}.log",
    params:
        locus=CHRY_LOCUS,
        vcf=CHRY_VCF,
    shell:
        """
        (
            bcftools view \
                --apply-filters PASS \
                --regions {params.locus} \
                --samples-file {input.samples_txt} \
                --force-samples \
                --output-type u \
                {params.vcf} \
            | bcftools annotate \
                --remove 'INFO,^FORMAT/GT' \
                --output-type z \
                --output {output.vcf} \
                --write-index=tbi \
                -
        ) &> {log}
        """


####################################################################################################
# Extract gnomAD allele frequencies for the chrY non-PAR test locus.
####################################################################################################
rule extract_gnomad_afs_chrY:
    input:
        variant_ht=f"{OUTPUT_DIR}/{CHRY_LOCUS_FILENAME}.ht",
    output:
        variant_ht=directory(f"{OUTPUT_DIR}/{CHRY_LOCUS_FILENAME}.gnomad_afs.ht"),
    log:
        f"logs/create_test_data/extract_gnomad_afs.{CHRY_LOCUS_FILENAME}.log",
    params:
        contig=CHRY_LOCUS_CHROM,
        freq_threshold=MIN_POP_AF_EXTRACT_GNOMAD_AFS,
        gcs_credentials_path=GCS_CREDENTIALS_PATH,
    shell:
        """
        (
            divref extract-gnomad-afs \
                --in-gnomad-sites-table {input.variant_ht} \
                --out-variant-annotation-table {output.variant_ht} \
                --contig {params.contig} \
                --freq-threshold {params.freq_threshold} \
                --gcs-credentials-path {params.gcs_credentials_path}
        ) &> {log}
        """
