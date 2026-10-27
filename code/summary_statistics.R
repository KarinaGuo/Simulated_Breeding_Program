## summary_statistics.R
##
## Calculates per-population summary statistics for each generation.
## Designed to be sourced and called from scenario scripts.
##
## Dependencies: tidyverse
##
## Usage:
##   source("summary_statistics.R")
##   stats <- calc_summary_stats(current_pop, original_pop)
##
## The `original_pop` (original_init_population) must contain:
##   $genotypes    - list of dosage vectors
##   $allele_freqs - ALT allele frequencies from SNP_AF.csv (the reference baseline)
##
## All diversity and capture statistics use original_pop$allele_freqs as the
## baseline, not the observed frequencies of the current population.

library(tidyverse)

## ── Helper: dosage matrix ─────────────────────────────────────────────────────
## Convert a population's genotype list to an N x L dosage matrix.
##   Rows = individuals, Columns = SNP loci.
pop_to_dosage_matrix <- function(pop) {
  do.call(rbind, lapply(pop$genotypes, as.integer))
}

## ── Allelic diversity ─────────────────────────────────────────────────────────
## Mean expected heterozygosity (He) across all loci.
##
## He_l = 2 * p_l * (1 - p_l)
##
## p_l is taken from original_pop$allele_freqs — the SNP_AF.csv baseline
## allele frequencies.  This ensures He is always measured against the same
## reference and is comparable across generations.
##
## He = 0   → all individuals fixed for the same allele
## He = 0.5 → maximum (allele frequency = 0.5)
##
## @param pop          Current population list
## @param allele_freqs Numeric vector of reference ALT allele freqs (from original_pop$allele_freqs)

calc_allelic_diversity <- function(pop, allele_freqs) {
  stopifnot(!is.null(allele_freqs), is.numeric(allele_freqs))
  He <- 2 * allele_freqs * (1 - allele_freqs)
  mean(He, na.rm = TRUE)
}

## ── Allelic capture ───────────────────────────────────────────────────────────
## Proportion of alleles present in `original_pop` that are also present in
## `current_pop` (i.e. dosage > 0 for ALT, or dosage < 2 for REF, in at least
## one individual of the current population).
##
## The SET of alleles to track is determined by original_pop$genotypes —
## alleles must actually appear in the F0 individuals to be counted.
## Allele presence/absence in current_pop is assessed from current_pop$genotypes.
##
## A value of 1 means no alleles have been lost relative to the F0.
##
## @param current_pop  Population list for the current generation
## @param original_pop The F0 population (original_init_population)

calc_allelic_capture <- function(current_pop, original_pop) {
  cur_dos  <- pop_to_dosage_matrix(current_pop)
  orig_dos <- pop_to_dosage_matrix(original_pop)

  n_loci   <- ncol(orig_dos)
  captured <- 0L
  total    <- 0L

  for (l in seq_len(n_loci)) {
    orig_col <- orig_dos[, l]
    cur_col  <- cur_dos[, l]

    ## REF allele present in F0?
    if (any(orig_col < 2)) {
      total    <- total + 1L
      if (any(cur_col < 2)) captured <- captured + 1L
    }
    ## ALT allele present in F0?
    if (any(orig_col > 0)) {
      total    <- total + 1L
      if (any(cur_col > 0)) captured <- captured + 1L
    }
  }

  if (total == 0L) return(NA_real_)
  captured / total
}

## ── Main summary function ─────────────────────────────────────────────────────
##
## @param current_pop   Population list for the current generation
## @param original_pop  The original F0 population (original_init_population).
##                      Must contain $allele_freqs for He computation.
##
## @return A named list with summary statistics

calc_summary_stats <- function(current_pop, original_pop) {

  stopifnot(
    is.list(current_pop),
    all(c("genotypes", "GEBV", "generation") %in% names(current_pop)),
    is.list(original_pop),
    all(c("genotypes", "GEBV", "allele_freqs") %in% names(original_pop))
  )

  ## Use original_pop$allele_freqs as the reference baseline for He
  list(
    generation        = current_pop$generation,
    n_individuals     = length(current_pop$indiv_ID),
    allelic_diversity = calc_allelic_diversity(current_pop,
                                               allele_freqs = original_pop$allele_freqs),
    allelic_capture   = calc_allelic_capture(current_pop, original_pop),
    mean_GEBV         = mean(current_pop$GEBV, na.rm = TRUE),
    sd_GEBV           = sd(current_pop$GEBV, na.rm = TRUE)
  )
}

## ── Pretty print helper ───────────────────────────────────────────────────────
print_summary_stats <- function(stats) {
  cat(sprintf(
    "── Summary Statistics (%s) ──────────────────────────\n",
    stats$generation
  ))
  cat(sprintf("  N individuals     : %d\n",   stats$n_individuals))
  cat(sprintf("  Allelic diversity : %.4f\n", stats$allelic_diversity))
  cat(sprintf("  Allelic capture   : %.4f\n", stats$allelic_capture))
  cat(sprintf("  Mean GEBV         : %.4f\n", stats$mean_GEBV))
  cat(sprintf("  SD GEBV           : %.4f\n", stats$sd_GEBV))
  cat("─────────────────────────────────────────────────────\n")
}

## ── Accumulate stats across generations ──────────────────────────────────────
## Usage:
##   stats_log <- list()
##   stats_log <- append_stats(stats_log, calc_summary_stats(f0_pop, f0_pop))
##   stats_log <- append_stats(stats_log, calc_summary_stats(f1_pop, f0_pop))
##   stats_df  <- stats_to_df(stats_log)

append_stats <- function(stats_log, stats) {
  c(stats_log, list(stats))
}

stats_to_df <- function(stats_log) {
  bind_rows(lapply(stats_log, as.data.frame))
}
