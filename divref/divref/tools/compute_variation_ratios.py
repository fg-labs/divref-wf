"""Tool to compute per-sample variant site counts across gnomAD frequency thresholds."""

import hail as hl

from divref import defaults
from divref.alias import HailPath


def compute_variation_ratios(
    *,
    vcfs_path: HailPath,
    gnomad_va_file: HailPath,
    gnomad_sa_file: HailPath,
    output_ht: HailPath,
    frequency_thresholds: list[float] = defaults.VARIATION_RATIO_FREQUENCY_THRESHOLDS,
) -> None:
    """
    Compute per-sample counts of non-reference variant sites across frequency thresholds.

    For each sample and each predefined gnomAD population frequency threshold, counts
    the number of sites where the sample carries a non-reference genotype. Writes a
    sample-level Hail table with population labels and per-threshold site counts.

    Args:
        vcfs_path: Path or glob pattern to input VCF files (GRCh38).
        gnomad_va_file: Path to the gnomAD variant annotations Hail table
            (from extract_gnomad_afs).
        gnomad_sa_file: Path to the gnomAD sample metadata Hail table
            (from extract_gnomad_afs).
        output_ht: Output path for the sample-level Hail table.
        frequency_thresholds: Frequency thresholds to calculate.
    """
    if any(t < 0 or t > 1 for t in frequency_thresholds):
        raise ValueError("All frequency_thresholds must be in [0, 1].")

    threshold_keys = [format(t, ".6g") for t in frequency_thresholds]
    if len(set(threshold_keys)) != len(threshold_keys):
        raise ValueError("frequency_thresholds produce duplicate output field names.")

    hl.init()

    gnomad_sa = hl.read_table(gnomad_sa_file)
    gnomad_va = hl.read_table(gnomad_va_file)

    mt = hl.import_vcf(vcfs_path, reference_genome="GRCh38", min_partitions=64, force_bgz=True)
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

    mt.cols().write(output_ht, overwrite=True)
