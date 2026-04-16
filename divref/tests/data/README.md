# Test data

The Hail table of gnomAD variant annotations at [chr1_100001_200000_0.01_seed42.ht](chr1_100001_200000_0.01_seed42.ht) contains a random sampling (seed 42) of 0.01 of the gnomAD HGDP+1KG variants in the chr1:100,001-200,000 locus from gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.hgdp_1kg_subset_variant_annotations.ht.

There are 91 variants in the table.

```bash
pixi run divref gnomad-hail-table-test-data \
  --out-variant-annotation-table chr1_100001_200000_0.01_seed42.ht
```
