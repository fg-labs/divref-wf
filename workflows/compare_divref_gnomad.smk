####################################################################################################
# Compares DivRef 1.1 single-variant sites against allele frequencies from different gnomAD
# releases for chr22.
####################################################################################################

from pathlib import Path

####################################################################################################
# Inputs / constants
####################################################################################################

OUTPUT_DIR: Path = Path("data/analysis")
CONTIG: str = "chr22"
DIVREF_DUCKDB_URL: str = (
    "https://zenodo.org/records/14802613/files/" "DivRef-v1.1.haplotypes_gnomad_merge.index.duckdb"
)
GNOMAD_VERSIONS: list[str] = ["joint_41", "hgdp_1kg_312"]

# Maps filename wildcard → plot label used by compare_divref_gnomad.R
GNOMAD_LABEL: dict[str, str] = {
    "joint_41": "gnomAD 4.1 joint AF",
    "hgdp_1kg_312": "gnomAD 3.1.2 HGDP+1KG AF",
}

####################################################################################################
# Rules
####################################################################################################


ruleorder: compare_divref_gnomad > extract_gnomad_single_afs


rule all:
    input:
        expand(
            f"{OUTPUT_DIR}/compare_divref_gnomad/{CONTIG}.{{gnomad_version}}.af_diffs.png",
            gnomad_version=GNOMAD_VERSIONS,
        ),
        expand(
            f"{OUTPUT_DIR}/compare_divref_gnomad/{CONTIG}.{{gnomad_version}}.not_in_gnomad_afs.png",
            gnomad_version=GNOMAD_VERSIONS,
        ),
        expand(
            f"{OUTPUT_DIR}/compare_divref_gnomad/{CONTIG}.{{gnomad_version}}.divref_not_in_gnomad.tsv",
            gnomad_version=GNOMAD_VERSIONS,
        ),
        expand(
            f"{OUTPUT_DIR}/compare_divref_gnomad/{CONTIG}.{{gnomad_version}}.log",
            gnomad_version=GNOMAD_VERSIONS,
        ),


####################################################################################################
# Downloads the DivRef 1.1 DuckDB index from Zenodo.
####################################################################################################
rule download_divref_index:
    output:
        duckdb=f"{OUTPUT_DIR}/input/DivRef-v1.1.haplotypes_gnomad_merge.index.duckdb",
    log:
        f"logs/compare_divref_gnomad/download_divref_index.log",
    params:
        url=DIVREF_DUCKDB_URL,
    shell:
        """
        (
            wget --no-verbose -O {output.duckdb} {params.url}
        ) &> {log}
        """


####################################################################################################
# Extracts chr22 allele frequencies from a gnomAD sites table and writes a Hail table and a flat
# TSV with one AF column per population.
####################################################################################################
rule extract_gnomad_single_afs:
    output:
        ht=directory(f"{OUTPUT_DIR}/compare_divref_gnomad/{CONTIG}.{{gnomad_version}}.ht"),
        tsv=f"{OUTPUT_DIR}/compare_divref_gnomad/{CONTIG}.{{gnomad_version}}.tsv",
    log:
        f"logs/compare_divref_gnomad/extract_gnomad_single_afs.{{gnomad_version}}.log",
    params:
        contig=CONTIG,
        gnomad_version=lambda wildcards: wildcards.gnomad_version.upper(),
    shell:
        """
        (
            divref extract-gnomad-single-afs \
                --contig {params.contig} \
                --gnomad-version {params.gnomad_version} \
                --out-sites-hail-table {output.ht} \
                --out-sites-tsv {output.tsv}
        ) &> {log}
        """


####################################################################################################
# Compares DivRef 1.1 gnomAD_variant sites to the extracted gnomAD allele frequencies.
# If all DivRef variants are found in gnomAD the R script exits early; touch ensures Snakemake's
# output checks are satisfied even in that case.
####################################################################################################
rule compare_divref_gnomad:
    input:
        duckdb=f"{OUTPUT_DIR}/input/DivRef-v1.1.haplotypes_gnomad_merge.index.duckdb",
        tsv=f"{OUTPUT_DIR}/compare_divref_gnomad/{CONTIG}.{{gnomad_version}}.tsv",
    output:
        af_diffs_png=f"{OUTPUT_DIR}/compare_divref_gnomad/{CONTIG}.{{gnomad_version}}.af_diffs.png",
        not_in_gnomad_png=f"{OUTPUT_DIR}/compare_divref_gnomad/{CONTIG}.{{gnomad_version}}.not_in_gnomad_afs.png",
        not_in_gnomad_tsv=f"{OUTPUT_DIR}/compare_divref_gnomad/{CONTIG}.{{gnomad_version}}.divref_not_in_gnomad.tsv",
    log:
        f"{OUTPUT_DIR}/compare_divref_gnomad/{CONTIG}.{{gnomad_version}}.log",
    params:
        contig=CONTIG,
        gnomad_label=lambda wildcards: GNOMAD_LABEL[wildcards.gnomad_version],
        output_base=f"{OUTPUT_DIR}/compare_divref_gnomad/{CONTIG}.{{gnomad_version}}",
    shell:
        """
        (
            Rscript scripts/compare_divref_gnomad.R \
                --contig {params.contig} \
                --divref_duckdb {input.duckdb} \
                --gnomad_tsv {input.tsv} \
                --gnomad_label '{params.gnomad_label}' \
                --output_base {params.output_base}

            # The R script exits early when no DivRef variants are absent from gnomAD,
            # so touch ensures these outputs exist for Snakemake's file checks.
            touch {output.not_in_gnomad_png} {output.not_in_gnomad_tsv}
        ) &> {log}
        """
