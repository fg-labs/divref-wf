"""Tool to export gnomAD variants above a frequency threshold as a VCF file."""

from pathlib import Path

import hail as hl

from divref.alias import HailPath
from divref.hail import hail_init


def create_gnomad_sites_vcf(
    *,
    sites_table_path: HailPath,
    output_vcf_path: HailPath,
    min_popmax: float,
    gcs_credentials_path: Path = Path("~/.config/gcloud/application_default_credentials.json"),
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
        gcs_credentials_path: Path to GCS default credentials JSON file.
    """
    hail_init(gcs_credentials_path.expanduser())

    ht = hl.read_table(sites_table_path)
    pops: list[str] = ht.pops.collect()[0]

    filt = ht.filter(hl.max(ht.pop_freqs.map(lambda x: x.AF)) >= min_popmax)
    filt = filt.annotate(
        info=hl.struct(**{
            f"{pop}_{fd}": filt.pop_freqs[i][fd]
            for i, pop in enumerate(pops)
            for fd in filt.pop_freqs[i]
        })
    )
    hl.export_vcf(filt, output_vcf_path)
