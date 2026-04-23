#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(duckdb)
  library(duckplyr)
  library(eulerr)
  library(logger)
  library(optparse)
  library(tidyverse)
})

option_list <- list(
  make_option("--contig",
    type = "character", default = "chr22",
    help = "Contig to filter from DivRef index [default: %default]"
  ),
  make_option("--divref_duckdb",
    type = "character",
    default = "data/resources/DivRef-v1.1.haplotypes_gnomad_merge.index.duckdb",
    help = "Path to DivRef DuckDB index [default: %default]"
  ),
  make_option("--gnomad_tsv",
    type = "character",
    help = "Path to gnomAD allele frequency TSV"
  ),
  make_option("--gnomad_label",
    type = "character",
    help = "Label for gnomAD data in plot axes and titles"
  ),
  make_option("--output_base",
    type = "character",
    help = "Base path for output files (suffixes added per file)"
  )
)

opts <- parse_args(OptionParser(option_list = option_list))

populations <- c("afr", "amr", "eas", "sas", "nfe")

# Load DivRef data ----

con <- dbConnect(duckdb(), read_only = TRUE, dbdir = opts$divref_duckdb)

divref <- dbGetQuery(
  con,
  "SELECT * FROM sequences WHERE contig = ? AND source = 'gnomAD_variant'",
  params = list(opts$contig)
) %>%
  select(-c(sequence, sequence_length, sequence_id, n_variants, source, contig)) %>%
  mutate(
    max_pop = as.factor(max_pop),
    gnomAD_AF_afr = as.numeric(gnomAD_AF_afr),
    gnomAD_AF_amr = as.numeric(gnomAD_AF_amr),
    gnomAD_AF_eas = as.numeric(gnomAD_AF_eas),
    gnomAD_AF_sas = as.numeric(gnomAD_AF_sas),
    gnomAD_AF_nfe = as.numeric(gnomAD_AF_nfe),
  ) %>%
  complete(fill = list(
    "gnomAD_AF_afr" = 0.0, "gnomAD_AF_amr" = 0.0, "gnomAD_AF_eas" = 0.0,
    "gnomAD_AF_sas" = 0.0, "gnomAD_AF_nfe" = 0.0
  ))

dbDisconnect(con)

log_info("Loaded ", nrow(divref), " DivRef variants for ", opts$contig)

# Load gnomAD data ----

gnomad <- read_tsv(opts$gnomad_tsv, show_col_types = FALSE) %>%
  mutate(
    afr = as.numeric(afr),
    amr = as.numeric(amr),
    eas = as.numeric(eas),
    sas = as.numeric(sas),
    nfe = as.numeric(nfe)
  ) %>%
  complete(fill = list("afr" = 0.0, "amr" = 0.0, "eas" = 0.0, "sas" = 0.0, "nfe" = 0.0))

# Join and split ----

divref_merged_with_gnomad <- divref %>%
  left_join(gnomad, by = join_by(variants == variant))

divref_in_gnomad <- divref_merged_with_gnomad %>%
  filter(!if_all(all_of(populations), is.na))

divref_not_in_gnomad <- divref_merged_with_gnomad %>%
  filter(if_all(all_of(populations), is.na)) %>%
  select(-all_of(populations))

log_info(nrow(divref_in_gnomad), " DivRef variants found in gnomAD")
log_info(nrow(divref_not_in_gnomad), " DivRef variants not found in gnomAD")

# Plot: Venn diagram of DivRef vs gnomAD variant overlap ----

n_both <- nrow(divref_in_gnomad)
n_divref_only <- nrow(divref_not_in_gnomad)
n_gnomad_only <- nrow(gnomad) - n_both

venn_counts <- c(n_divref_only, n_gnomad_only, n_both)
names(venn_counts) <- c("DivRef 1.1", opts$gnomad_label, paste0("DivRef 1.1&", opts$gnomad_label))
fit <- euler(venn_counts)

png(paste0(opts$output_base, ".venn.png"), height = 600, width = 600)
g <- plot(fit, quantities = TRUE)
grid::grid.newpage()
grid::pushViewport(grid::viewport(width = 0.8, height = 0.8))
grid::grid.draw(g)
grid::popViewport()
dev.off()

# Plot: AF differences for variants found in gnomAD ----

divref_in_gnomad_with_af_diffs <- divref_in_gnomad %>%
  mutate(
    diff_afr = afr - gnomAD_AF_afr,
    diff_amr = amr - gnomAD_AF_amr,
    diff_eas = eas - gnomAD_AF_eas,
    diff_sas = sas - gnomAD_AF_sas,
    diff_nfe = nfe - gnomAD_AF_nfe,
  )

p <- divref_in_gnomad_with_af_diffs %>%
  select(variants, diff_afr, diff_amr, diff_eas, diff_sas, diff_nfe) %>%
  pivot_longer(
    cols = c(diff_afr, diff_amr, diff_eas, diff_sas, diff_nfe),
    names_to = "population", values_to = "diff_freq", names_prefix = "diff_"
  ) %>%
  ggplot(aes(x = diff_freq)) +
  geom_histogram() +
  facet_wrap(~population, nrow = 5) +
  scale_y_log10() +
  theme_bw() +
  xlab(paste0(opts$gnomad_label, " - DivRef 1.1 AF")) +
  ylab("Variants")

ggsave(paste0(opts$output_base, ".af_diffs.png"), p, height = 10, width = 6)

p <- divref_in_gnomad_with_af_diffs %>%
  select(variants, diff_afr, diff_amr, diff_eas, diff_sas, diff_nfe) %>%
  pivot_longer(
    cols = c(diff_afr, diff_amr, diff_eas, diff_sas, diff_nfe),
    names_to = "population", values_to = "diff_freq", names_prefix = "diff_"
  ) %>%
  ggplot(aes(x = diff_freq)) +
  geom_histogram() +
  scale_y_log10() +
  theme_bw() +
  xlab(paste0(opts$gnomad_label, " - DivRef 1.1 AF")) +
  ylab("Variants")

ggsave(paste0(opts$output_base, ".af_diffs_all.png"), p, height = 6, width = 6)

# Count: variants with large AF differences ----

n_large_af_diff <- divref_in_gnomad_with_af_diffs %>%
  select(variants, diff_afr, diff_amr, diff_eas, diff_sas, diff_nfe) %>%
  pivot_longer(
    cols = c(diff_afr, diff_amr, diff_eas, diff_sas, diff_nfe),
    names_to = "population", values_to = "diff_freq", names_prefix = "diff_"
  ) %>%
  mutate(abs_diff_freq = abs(diff_freq)) %>%
  filter(abs_diff_freq >= 0.001) %>%
  distinct(variants) %>%
  nrow()

log_info(n_large_af_diff, " DivRef variants found in gnomAD with |AF diff| >= 0.001 in any population")

if (nrow(divref_not_in_gnomad) == 0) {
  log_info("No DivRef variants missing from gnomAD, exiting")
  quit(save = "no")
}

# Plot: DivRef AF distribution for variants not found in gnomAD ----

p <- divref_not_in_gnomad %>%
  select(variants, gnomAD_AF_afr, gnomAD_AF_amr, gnomAD_AF_eas, gnomAD_AF_sas, gnomAD_AF_nfe) %>%
  pivot_longer(
    cols = c(gnomAD_AF_afr, gnomAD_AF_amr, gnomAD_AF_eas, gnomAD_AF_sas, gnomAD_AF_nfe),
    names_to = "population", values_to = "freq", names_prefix = "gnomAD_AF_"
  ) %>%
  ggplot(aes(x = freq)) +
  geom_histogram() +
  facet_wrap(~population, nrow = 5) +
  scale_y_log10() +
  theme_bw() +
  xlab("DivRef 1.1 AF") +
  ylab("Variants") +
  labs(title = paste0(
    "DivRef 1.1 'gnomAD_variant' variants not found in ", opts$gnomad_label
  ))

ggsave(paste0(opts$output_base, ".not_in_gnomad_afs.png"), p, height = 10, width = 6)

# Write not-in-gnomAD variant list ----

divref_not_in_gnomad %>%
  select(variants) %>%
  write_tsv(paste0(opts$output_base, ".divref_not_in_gnomad.tsv"))
