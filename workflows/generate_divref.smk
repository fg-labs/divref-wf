####################################################################################################
# Generates a DivRef-format resource of human haplotypes.
#
# Final output is a set of per-chromosome FASTA files and a DuckDB index.
####################################################################################################

from pathlib import Path
from snakemake.utils import validate

####################################################################################################
# Inputs
####################################################################################################


configfile: os.path.join(workflow.basedir, "config", "config.yml")


validate(config, os.path.join(workflow.basedir, "config", "config_schema.yml"))

WORK_DIR: Path = Path(config["work_dir"])
CHROMS: list[str] = config["chromosomes"]

# The HGDP+1KG phased BCF files are at
# "{HGDP_1KG_PHASED_BCF_PREFIX}.{chrom}.{HGDP_1KG_PHASED_BCF_SUFFIX}"
HGDP_1KG_PHASED_BCF_PREFIX: str = config["hgdp_1kg_phased_bcf_prefix"]
HGDP_1KG_PHASED_BCF_SUFFIX: str = config["hgdp_1kg_phased_bcf_suffix"]

VCF_EXTS: list[str] = [".vcf.gz", ".vcf.gz.tbi"]

####################################################################################################
# Rules
####################################################################################################


rule all:
    input:
        #        f"{WORK_DIR}/inputs/hgdp_1kg_sample_metadata.ht",
        #expand(f"{WORK_DIR}/inputs/hgdp_1kg.sites.{{chrom}}.ht", chrom=CHROMS),
        expand(
            f"{WORK_DIR}/inputs/hgdp_1kg.phased_genotypes.{{chrom}}{{ext}}",
            chrom=CHROMS,
            ext=VCF_EXTS,
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
