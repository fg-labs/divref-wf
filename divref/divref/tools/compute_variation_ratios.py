"""Tool to compute per-sample variant site counts across gnomAD frequency thresholds."""

import hail as hl

from divref.alias import HailPath

_FREQ_THRESHOLDS = [0, 0.0001, 0.001, 0.005, 0.01, 0.05, 0.1]


def compute_variation_ratios(
    *,
    vcfs_path: HailPath,
    gnomad_va_file: HailPath,
    gnomad_sa_file: HailPath,
    output_ht: HailPath,
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
    """
    hl.init()

    gnomad_sa = hl.read_table(gnomad_sa_file)
    gnomad_va = hl.read_table(gnomad_va_file)

    mt = hl.import_vcf(vcfs_path, reference_genome="GRCh38", min_partitions=64)
    mt = mt.select_rows().select_cols()
    mt = mt.annotate_rows(freq=gnomad_va[mt.row_key].pop_freqs)
    mt = mt.filter_rows(hl.is_defined(mt.freq))
    mt = mt.annotate_rows(max_pop_freq=hl.max(mt.freq.map(lambda x: hl.max(x.AF))))

    mt = mt.annotate_cols(pop=gnomad_sa[mt.col_key].pop)
    mt = mt.annotate_cols(
        counts=hl.struct(**{
            f"n_sites_above_{x}": hl.agg.count_where(mt.GT.is_non_ref() & (mt.max_pop_freq > x))
            for x in _FREQ_THRESHOLDS
        })
    )

    mt.cols().write(output_ht, overwrite=True)
