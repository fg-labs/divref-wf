"""Tool to compute haplotype and gnomAD variant statistics across parameter combinations."""

import hail as hl
from pydantic import BaseModel

from divref.haplotype import HailPath
from divref.haplotype import split_haplotypes


class _HGDPResult(BaseModel):
    """Summary of unique haplotype count for one (frequency, window_size) parameter combination."""

    frequency_cutoff: float
    window_size: int
    hgdp_haplotype_count: int


class _GnomADResult(BaseModel):
    """Count of gnomAD variants above a given frequency cutoff."""

    frequency_cutoff: float
    gnomad_variant_count: int


def compute_haplotype_statistics(
    *,
    haplotypes_table_path: HailPath,
    gnomad_va_file: HailPath,
    window_sizes: list[int],
    frequency_cutoffs: list[float],
    output_base: HailPath,
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
    """
    hl.init()

    ht = hl.read_table(haplotypes_table_path).key_by()
    va = hl.read_table(gnomad_va_file)

    hgdp_results: list[_HGDPResult] = []
    gnomad_results: list[_GnomADResult] = []

    for frequency in frequency_cutoffs:
        for window_size in window_sizes:
            ht2 = ht.filter(ht.max_pop_freq >= frequency)
            ht2 = split_haplotypes(ht2, window_size)
            n_unique = ht2.key_by("haplotype").distinct().key_by().count()
            hgdp_results.append(
                _HGDPResult(
                    frequency_cutoff=frequency,
                    window_size=window_size,
                    hgdp_haplotype_count=n_unique,
                )
            )

        gnomad_count = va.filter(
            hl.max(va.pop_freqs.map(lambda x: hl.max(x.AF))) >= frequency
        ).count()
        gnomad_results.append(
            _GnomADResult(frequency_cutoff=frequency, gnomad_variant_count=gnomad_count)
        )

    with open(f"{output_base}.hgdp.tsv", "w") as f:
        f.write("frequency\twindow_size\thgdp_haplotype_count\n")
        for hgdp_result in hgdp_results:
            f.write(
                f"{hgdp_result.frequency_cutoff}\t"
                f"{hgdp_result.window_size}\t"
                f"{hgdp_result.hgdp_haplotype_count}\n"
            )

    with open(f"{output_base}.gnomad.tsv", "w") as f:
        f.write("frequency\tgnomad_variant_count\n")
        for gnomad_result in gnomad_results:
            f.write(f"{gnomad_result.frequency_cutoff}\t{gnomad_result.gnomad_variant_count}\n")
