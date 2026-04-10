from snakemake.utils import validate

configfile: os.path.join(workflow.basedir, "config", "config.yml")

validate(config, os.path.join(workflow.basedir, "config", "config_schema.yml"))

SAMPLES = config["samples"]
READS = ["1", "2"]


rule all:
    input:
        "data/resources/ref.fa",
        expand("data/raw/{sample}_R{read}.fastq.gz", sample=SAMPLES, read=READS),


rule download_raw_data:
    output:
        "data/raw/{sample}_R{read}.fastq.gz",
    shell:
        """
        # wget -O {output} https://to/data/{wildcards.sample}_R{wildcards.read}.fastq.gz
        touch {output}
        """


rule download_resource_data:
    output:
        "data/resources/ref.fa",
    shell:
        """
        # wget -O {output} https://to/data/ref.fa
        touch {output}
        """
