"""Tool to compute per-sample variant site counts across gnomAD frequency thresholds."""

import os
from pathlib import Path

import hail as hl
from fgpyo.io import assert_directory_exists
from fgpyo.io import assert_path_is_readable
from fgpyo.io import assert_path_is_writable

from divref import defaults


def compute_variation_ratios(
    *,
    vcfs_path: Path,
    gnomad_va_file: Path,
    gnomad_sa_file: Path,
    output_ht: Path,
    frequency_thresholds: list[float] = defaults.VARIATION_RATIO_FREQUENCY_THRESHOLDS,
    reference_genome: str = defaults.REFERENCE_GENOME,
    spark_driver_memory_gb: int = 1,
    spark_executor_memory_gb: int = 1,
) -> None:
    """
    Compute per-sample counts of non-reference variant sites across frequency thresholds.

    For each sample and each predefined gnomAD population frequency threshold, counts
    the number of sites where the sample carries a non-reference genotype. Writes a
    sample-level Hail table with population labels and per-threshold site counts.

    Args:
        vcfs_path: Path to input VCF files.
        gnomad_va_file: Path to the gnomAD variant annotations Hail table
            (from extract_gnomad_afs).
        gnomad_sa_file: Path to the gnomAD sample metadata Hail table
            (from extract_gnomad_afs).
        output_ht: Output path for the sample-level Hail table.
        frequency_thresholds: Frequency thresholds to calculate.
        reference_genome: Reference genome to use. Defaults to "GRCh38".
        spark_driver_memory_gb: Memory in GB to allocate to the Spark driver.
        spark_executor_memory_gb: Memory in GB to allocate to the Spark executor.
    """
    assert_path_is_readable(vcfs_path)
    assert_directory_exists(gnomad_va_file)
    assert_directory_exists(gnomad_sa_file)
    assert_path_is_writable(output_ht)

    if any(t < 0 or t > 1 for t in frequency_thresholds):
        raise ValueError("All frequency_thresholds must be in [0, 1].")

    threshold_keys = [format(t, ".6g") for t in frequency_thresholds]
    if len(set(threshold_keys)) != len(threshold_keys):
        raise ValueError("frequency_thresholds produce duplicate output field names.")

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

    gnomad_sa = hl.read_table(str(gnomad_sa_file))
    gnomad_va = hl.read_table(str(gnomad_va_file))

    mt = hl.import_vcf(str(vcfs_path), reference_genome=reference_genome, min_partitions=64)
    mt = mt.select_rows().select_cols()
    mt = mt.annotate_rows(freq=gnomad_va[mt.row_key].pop_freqs)
    mt = mt.filter_rows(hl.is_defined(mt.freq))
    mt = mt.annotate_rows(max_pop_freq=hl.max(mt.freq.map(lambda x: x.AF)))

    mt = mt.annotate_cols(pop=gnomad_sa[mt.col_key].pop)
    mt = mt.annotate_cols(
        counts=hl.struct(**{
            f"n_sites_above_{key.replace('.', '_')}": hl.agg.count_where(
                mt.GT.is_non_ref() & (mt.max_pop_freq > threshold)
            )
            for threshold, key in zip(frequency_thresholds, threshold_keys, strict=True)
        })
    )

    mt.cols().write(str(output_ht), overwrite=True)
