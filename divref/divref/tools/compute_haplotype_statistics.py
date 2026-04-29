"""Tool to compute haplotype and gnomAD variant statistics across parameter combinations."""

import os
from pathlib import Path

import hail as hl
from fgmetric import Metric
from fgmetric import MetricWriter
from fgpyo.io import assert_directory_exists
from fgpyo.io import assert_path_is_writable

from divref.haplotype import split_haplotypes


class _HGDPResult(Metric):
    """Summary of unique haplotype count for one (frequency, window_size) parameter combination."""

    frequency_cutoff: float
    window_size: int
    hgdp_haplotype_count: int


class _GnomADResult(Metric):
    """Count of gnomAD variants above a given frequency cutoff."""

    frequency_cutoff: float
    gnomad_variant_count: int


def compute_haplotype_statistics(
    *,
    haplotypes_table_path: Path,
    gnomad_va_file: Path,
    window_sizes: list[int],
    frequency_cutoffs: list[float],
    output_base: Path,
    spark_driver_memory_gb: int = 1,
    spark_executor_memory_gb: int = 1,
) -> None:
    """
    Compute haplotype and gnomAD variant statistics across frequency and window size parameters.

    For each combination of frequency cutoff and window size, counts unique haplotypes in
    the HGDP dataset. For each frequency cutoff, counts gnomAD variants above that threshold.
    Writes two TSV summary files: one for haplotype counts and one for gnomAD variant counts.

    Args:
        haplotypes_table_path: Path to the Hail table of computed haplotypes
            (from compute_haplotypes).
        gnomad_va_file: Path to the gnomAD variant annotations Hail table
            (from extract_gnomad_afs).
        window_sizes: Window sizes in bp to evaluate (e.g. --window-sizes 100 500 1000).
        frequency_cutoffs: Allele frequency thresholds to evaluate
            (e.g. --frequency-cutoffs 0.001 0.005 0.01).
        output_base: Base path for output TSV files; writes {output_base}.hgdp.tsv
            and {output_base}.gnomad.tsv.
        spark_driver_memory_gb: Memory in GB to allocate to the Spark driver.
        spark_executor_memory_gb: Memory in GB to allocate to the Spark executor.
    """
    assert_directory_exists(haplotypes_table_path)
    assert_directory_exists(gnomad_va_file)

    out_hgdp: Path = output_base.with_suffix(".hgdp.tsv")
    out_gnomad: Path = output_base.with_suffix(".gnomad.tsv")
    assert_path_is_writable(out_hgdp)
    assert_path_is_writable(out_gnomad)

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

    ht = hl.read_table(str(haplotypes_table_path)).key_by()
    va = hl.read_table(str(gnomad_va_file))

    with (
        MetricWriter(_HGDPResult, out_hgdp) as hgdp_writer,
        MetricWriter(_GnomADResult, out_gnomad) as gnomad_writer,
    ):
        for frequency in frequency_cutoffs:
            ht_filtered = ht.filter(
                hl.max(ht.all_pop_freqs.map(lambda x: x.empirical_AF)) >= frequency
            )
            for window_size in window_sizes:
                ht2 = split_haplotypes(ht_filtered, window_size)
                n_unique = ht2.key_by("haplotype").distinct().key_by().count()
                hgdp_writer.write(
                    _HGDPResult(
                        frequency_cutoff=frequency,
                        window_size=window_size,
                        hgdp_haplotype_count=n_unique,
                    )
                )

            gnomad_count = va.filter(hl.max(va.pop_freqs.map(lambda x: x.AF)) >= frequency).count()
            gnomad_writer.write(
                _GnomADResult(frequency_cutoff=frequency, gnomad_variant_count=gnomad_count)
            )
