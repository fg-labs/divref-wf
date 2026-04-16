"""Tool to compute haplotypes from VCF files with gnomAD population frequency annotations."""

import logging
from typing import Any
from typing import Callable

import hail as hl

from divref.alias import HailPath

logger = logging.getLogger(__name__)


def _get_haplotypes(
    ht: Any,
    windower_f: Callable[[Any], Any],
    idx: int,
    output_base: HailPath,
    pop_ints: dict[str, int],
) -> Any:
    """
    Group variants into haplotypes within genomic windows and compute empirical frequencies.

    Applies the windowing function to assign variants to windows, aggregates haplotypes per
    population and sample, filters to multi-variant haplotypes, and collapses across samples
    to compute empirical allele counts and frequencies. Writes intermediate results to a
    checkpoint table.

    Args:
        ht: Hail table with per-variant population membership and frequency data.
        windower_f: Function mapping a Hail locus to the window locus key.
        idx: Index of this windowing pass (1 or 2), used in the checkpoint filename.
        output_base: Base path for output; checkpoint written to {output_base}.{idx}.ht.
        pop_ints: Mapping from population code to integer index.

    Returns:
        Hail table of haplotypes with empirical frequency summaries.
    """
    new_locus = windower_f(ht.locus)
    ht = ht.annotate(new_locus=new_locus)

    def agg_haplos(arr: Any) -> Any:
        flat = hl.agg.explode(lambda elt: hl.agg.collect(elt.annotate(row_idx=ht.row_idx)), arr)
        pop_grouped = hl.group_by(lambda x: x.pop, flat)
        return pop_grouped.map_values(
            lambda arr_per_pop: hl.array(
                hl.array(hl.group_by(lambda inner_elt: inner_elt.sample, arr_per_pop))
                .filter(lambda sample_and_records: hl.len(sample_and_records[1]) > 1)
                .map(
                    lambda sample_and_records: hl.sorted(
                        sample_and_records[1].map(lambda e: e.row_idx)
                    )
                )
                .group_by(lambda x: x)
                .map_values(lambda arr: hl.len(arr))
            )
        )

    ht_grouped = ht.group_by("new_locus").aggregate(
        row_map=hl.dict(
            hl.agg.collect((
                ht.row_idx,
                ht.row.select("locus", "alleles", "freq", "frequencies_by_pop"),
            ))
        ),
        left_haplos=agg_haplos(ht.pops_and_ids_left),
        right_haplos=agg_haplos(ht.pops_and_ids_right),
    )

    def collapse_haplos_across_samples(pop: Any, arr1: Any, arr2: Any) -> Any:
        # Assumes all AN == 2 * N_samples.
        flat = hl.array([arr1, arr2]).flatmap(lambda x: x.get(pop))

        def map_haplo_group(t: Any) -> Any:
            haplotype = t[0]
            n_observed = hl.sum(t[1].map(lambda x: x[1]))
            component_variant_frequencies = haplotype.map(
                lambda x: ht_grouped.row_map[x].frequencies_by_pop[pop]
            )
            min_an = hl.min(component_variant_frequencies.map(lambda x: x.AN))
            return hl.struct(
                haplotype=haplotype,
                pop=pop,
                empirical_AC=n_observed,
                min_variant_frequency=hl.min(component_variant_frequencies.map(lambda x: x.AF[1])),
                empirical_AF=n_observed / min_an,
            )

        return hl.array(hl.group_by(lambda x: x[0], flat)).map(map_haplo_group)

    ht_grouped = ht_grouped.annotate(
        all_haplos=hl.literal(list(pop_ints.values())).flatmap(
            lambda pop: collapse_haplos_across_samples(
                pop, ht_grouped.left_haplos, ht_grouped.right_haplos
            )
        )
    )

    def get_haplotype_summary(a: Any) -> dict[str, Any]:
        a_sorted = hl.sorted(a, key=lambda x: x.empirical_AF, reverse=True)
        return dict(
            max_pop=a_sorted[0].pop,
            max_empirical_AF=a_sorted[0].empirical_AF,
            max_empirical_AC=a_sorted[0].empirical_AC,
            min_variant_frequency=a_sorted[0].min_variant_frequency,
            all_pop_freqs=a_sorted.map(lambda x: x.drop("haplotype")),
        )

    ht_grouped = ht_grouped.transmute(
        all_haplos=hl.array(hl.group_by(lambda x: x.haplotype, ht_grouped.all_haplos)).map(
            lambda t: hl.struct(haplotype=t[0], **get_haplotype_summary(t[1]))
        )
    )

    hte = ht_grouped.explode("all_haplos")
    hte = hte.key_by().drop("new_locus")

    def get_variant(row_idx: Any) -> Any:
        return hte.row_map[row_idx].select("locus", "alleles")

    def get_gnomad_freq(row_idx: Any) -> Any:
        return hte.row_map[row_idx].freq

    hte = hte.select(
        **hte.all_haplos,
        variants=hte.all_haplos.haplotype.map(get_variant),
        gnomad_freqs=hte.all_haplos.haplotype.map(get_gnomad_freq),
    )

    hte = hte.group_by("haplotype").aggregate(
        **hl.sorted(
            hl.agg.collect(hte.row.drop("haplotype")),
            key=lambda row: -row.max_empirical_AF,
        )[0]
    )

    logger.info("Writing %s.%s.ht ...", output_base, idx)
    return hte.checkpoint(f"{output_base}.{idx}.ht", overwrite=True)


def compute_haplotypes(
    *,
    vcfs_path: HailPath,
    gnomad_va_file: HailPath,
    gnomad_sa_file: HailPath,
    window_size: int,
    freq_threshold: float,
    output_base: HailPath,
    temp_dir: HailPath = "/tmp",
) -> None:
    """
    Compute population haplotypes from VCF files with gnomAD frequency annotations.

    Reads VCF files, annotates variants with gnomAD population allele frequencies,
    extracts phased haplotypes per population using two overlapping window strategies,
    and writes the union of both windowed results as a keyed Hail table.

    Args:
        vcfs_path: Path or glob pattern to input VCF files (GRCh38).
        gnomad_va_file: Path to the gnomAD variant annotations Hail table
            (from extract_gnomad_afs).
        gnomad_sa_file: Path to the gnomAD sample metadata Hail table
            (from extract_gnomad_afs).
        window_size: Base window size in bp for grouping variants into haplotypes.
        freq_threshold: Minimum gnomAD population allele frequency to retain a variant.
        output_base: Base output path; writes {output_base}.1.ht, {output_base}.2.ht,
            and the final {output_base}.ht.
        temp_dir: Directory for Hail temporary files.
    """
    hl.init(tmp_dir=temp_dir)

    gnomad_sa = hl.read_table(gnomad_sa_file)
    gnomad_va = hl.read_table(gnomad_va_file)
    gnomad_va = gnomad_va.filter(hl.max(gnomad_va.pop_freqs.map(lambda x: x.AF)) >= freq_threshold)

    mt = hl.import_vcf(vcfs_path, reference_genome="GRCh38", min_partitions=64)
    mt = mt.select_rows().select_cols()
    mt = mt.annotate_rows(freq=gnomad_va[mt.row_key].pop_freqs)
    mt = mt.filter_rows(hl.is_defined(mt.freq))

    pop_legend: list[str] = gnomad_va.globals.pops.collect()[0]
    pop_ints = {pop: i for i, pop in enumerate(pop_legend)}
    mt = mt.annotate_cols(pop_int=hl.literal(pop_ints).get(gnomad_sa[mt.col_key].pop))
    mt = mt.filter_cols(hl.is_defined(mt.pop_int))
    mt = mt.add_row_index().add_col_index()
    mt = mt.filter_entries(mt.freq[mt.pop_int].AF >= freq_threshold)

    mt = mt.annotate_rows(
        pops_and_ids_left=hl.agg.filter(
            mt.GT[0] != 0, hl.agg.collect(hl.struct(pop=mt.pop_int, sample=mt.col_idx))
        ),
        pops_and_ids_right=hl.agg.filter(
            mt.GT[1] != 0, hl.agg.collect(hl.struct(pop=mt.pop_int, sample=mt.col_idx))
        ),
        frequencies_by_pop=hl.agg.group_by(mt.pop_int, hl.agg.call_stats(mt.GT, 2)),
    )
    ht = mt.rows().select(
        "freq",
        "pops_and_ids_left",
        "pops_and_ids_right",
        "row_idx",
        "frequencies_by_pop",
    )

    window1 = _get_haplotypes(
        ht,
        lambda locus: locus - (locus.position % window_size),
        1,
        output_base,
        pop_ints,
    )
    window2 = _get_haplotypes(
        ht,
        lambda locus: locus - ((locus.position + window_size // 2) % window_size),
        2,
        output_base,
        pop_ints,
    )

    htu = window1.union(window2)
    logger.info("Writing final %s.ht ...", output_base)
    htu.key_by("haplotype").naive_coalesce(64).write(f"{output_base}.ht", overwrite=True)
