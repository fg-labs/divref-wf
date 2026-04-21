####################################################################################################
# Creates the test data in divref/tests/data.
####################################################################################################

from pathlib import Path

####################################################################################################
# Inputs
####################################################################################################

OUTPUT_DIR: Path = Path("data/test")
LOCUS_CHROM: str = "chr1"
LOCUS: str = "chr1:100001-200000"
LOCUS_FILENAME: str = "chr1_100001_200000"

####################################################################################################
# Rules
####################################################################################################


rule all:
    input:
        f"{OUTPUT_DIR}/{LOCUS_FILENAME}.ht",
        f"{OUTPUT_DIR}/hgdp_1kg_sample_metadata.ht",
        f"{OUTPUT_DIR}/{LOCUS_FILENAME}.vcf.gz",
        f"{OUTPUT_DIR}/{LOCUS_FILENAME}.vcf.gz.tbi",


####################################################################################################
# Extracts all gnomAD HGDP+1KG variants in the specified locus from
# gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_variant_annotations.ht.
#
# Extracts selected fields from gnomAD sample metadata required by downstream tools from
# gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.ht.
####################################################################################################
rule subset_gnomad_hail_tables:
    output:
        variant_ht=directory(f"{OUTPUT_DIR}/{LOCUS_FILENAME}.ht"),
        sample_ht=directory(f"{OUTPUT_DIR}/hgdp_1kg_sample_metadata.ht"),
    log:
        f"logs/create_test_data/subset_gnomad_hail_tables.{LOCUS_FILENAME}.log",
    params:
        locus=LOCUS,
    shell:
        """
        (
            divref gnomad-hail-table-test-data \
                --out-variant-annotation-table {output.variant_ht} \
                --out-sample-metadata {output.sample_ht} \
                --locus {params.locus}
        ) &> {log}
        """


####################################################################################################
# Extracts the phased genotypes for all HGDP+1KG samples in the specified locus.
####################################################################################################
rule subset_phased_genotypes:
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
                --output-type z \
                --output {output.vcf} \
                --write-index=tbi \
                {params.bcf}
        ) &> {log}
        """
