## scenario_utils.R
## Shared utilities and setup for all project scenarios.
setwd("~/Uni/Doctorate/Breeding_project")
library(tidyverse)

## Source dependencies
source("code/gapit_functions_080425.txt")
source("code/GP_functs.R")
source("code/genotype_phenotype_v1.R")
source("code/summary_statistics.R")
source("code/calculate_kinship.R")
source("code/optimise_crosses.R")
source("code/optimise_grid_placement.R")

## Setup global variables required by GAPIT pipeline
out_dir             <<- "output/"
gt_datafile         <<- "data/Report-DMela25-10229/Report_DMela25-10229_RegularGenotyping/Report_DMela25-10229_GenotypingSamples_trainingconcat_sort.hapmap.hmp.txt"
input_phenodatafile <<- "data/mq_phenotypes.csv"

load("data/original_data.rds")

## Load SNP lists
#SNP_list    <<- read.csv("data/SNP_list.csv", header = TRUE)
#SNPs_tested <<- SNP_list[, 1]

## Load metadata from training hapmap
df_base <- read.csv(gt_datafile, sep = "\t", header = TRUE) %>%
  filter(rs. %in% SNPs_tested)

hapmap_meta <<- df_base[, 1:11]
training_hapmap <<- df_base[, c(1, grep("S_", colnames(df_base)))]

## Establish SNP REF/ALT and IUPAC codes
snp_refalt_identity <- df_base %>%
  separate(col = alleles, into = c("REF", "ALT"), sep = "/", fill = "right", remove = FALSE) %>%
  dplyr::select(rs., REF, ALT)

SNPs_tested_ord <- match(SNPs_tested, snp_refalt_identity$rs.)
snp_refalt_identity_reord <- snp_refalt_identity[SNPs_tested_ord, ]
ref_alleles <<- snp_refalt_identity_reord$REF
alt_alleles <<- snp_refalt_identity_reord$ALT

iupac_dict <- c("AG" = "R", "GA" = "R", "CT" = "Y", "TC" = "Y", "GC" = "S", "CG" = "S", 
                "AT" = "W", "TA" = "W", "GT" = "K", "TG" = "K", "AC" = "M", "CA" = "M")
het_codes <<- iupac_dict[paste0(ref_alleles, alt_alleles)]

## Load Initial Population
original_init_population <<- readRDS("data/original_init_population.rds")

## Helper: Subset a population to the top fraction by GEBV
get_top_fraction <- function(pop, fraction) {
  n_keep <- max(1, floor(length(pop$indiv_ID) * fraction))
  idx <- order(pop$GEBV)[1:n_keep]
  list(
    indiv_ID    = pop$indiv_ID[idx],
    genotypes   = pop$genotypes[idx],
    GEBV        = pop$GEBV[idx],
    allele_freqs = pop$allele_freqs,
    maternal_ID = pop$maternal_ID[idx],
    paternal_ID = pop$paternal_ID[idx],
    generation  = pop$generation,
    error       = pop$error[idx]
  )
}

## Helper: Generate a random crossing plan (repeats allowed)
random_crossing_plan <- function(maternal_pop, paternal_pop, target_n = 100, allow_selfing = TRUE) {
  mat_ids <- as.character(maternal_pop$indiv_ID)
  pat_ids <- as.character(paternal_pop$indiv_ID)
  
  n_mat <- length(mat_ids)
  n_each <- floor(target_n / n_mat)
  n_extra <- target_n - (n_each * n_mat)
  
  plan <- vector("list", n_mat)
  for (i in seq_len(n_mat)) {
    candidates <- if (!allow_selfing) setdiff(pat_ids, mat_ids[i]) else pat_ids
    pat_id <- sample(candidates, 1)
    
    n_off <- n_each + if (i <= n_extra) 1L else 0L
    plan[[i]] <- data.frame(
      maternal_ID = mat_ids[i],
      paternal_ID = pat_id,
      n_offspring = n_off,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, plan)
}
