"""Tool to export gnomAD variants above a frequency threshold as a VCF file."""

import os
from pathlib import Path

import hail as hl
from fgpyo.io import assert_directory_exists
from fgpyo.io import assert_path_is_writable


def create_gnomad_sites_vcf(
    *,
    sites_table_path: Path,
    output_vcf_path: Path,
    min_popmax: float,
    spark_driver_memory_gb: int = 1,
    spark_executor_memory_gb: int = 1,
) -> None:
    """
    Export gnomAD variant sites above a population frequency threshold as VCF.

    Filters the gnomAD variant annotations table to sites where at least one
    population's allele frequency meets or exceeds min_popmax, restructures
    per-population frequency fields into VCF INFO format, and exports as VCF.

    Args:
        sites_table_path: Path to the gnomAD variant annotations Hail table output from
            `extract-gnomad-afs`.
        output_vcf_path: Output path for the VCF file.
        min_popmax: Minimum allele frequency in any population to include a site.
        spark_driver_memory_gb: Memory in GB to allocate to the Spark driver.
        spark_executor_memory_gb: Memory in GB to allocate to the Spark executor.
    """
    assert_directory_exists(sites_table_path)
    assert_path_is_writable(output_vcf_path)

    if spark_driver_memory_gb < 1:
        raise ValueError(
            f"Spark driver memory must be at least 1GB. Saw {spark_driver_memory_gb}GB."
        )
    if spark_executor_memory_gb < 1:
        raise ValueError(
            f"Spark executor memory must be at least 1GB. Saw {spark_executor_memory_gb}GB."
        )

    os.environ["PYSPARK_SUBMIT_ARGS"] = (
        f"--driver-memory {spark_driver_memory_gb}g "
        f"--executor-memory {spark_executor_memory_gb}g "
        "pyspark-shell"
    )
    hl.init()

    ht = hl.read_table(str(sites_table_path))
    pops: list[str] = ht.pops.collect()[0]

    filt = ht.filter(hl.max(ht.pop_freqs.map(lambda x: x.AF)) >= min_popmax)
    filt = filt.annotate(
        info=hl.struct(**{
            f"{pop}_{fd}": filt.pop_freqs[i][fd]
            for i, pop in enumerate(pops)
            for fd in filt.pop_freqs[i]
        })
    )
    hl.export_vcf(filt, str(output_vcf_path))
