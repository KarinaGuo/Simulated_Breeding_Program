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

## Load Fixed Maternal Individuals if provided (optional)
## Expected format: A dataframe/CSV with columns indiv_ID, and SNP dosages (0,1,2) for all SNPs_tested
fixed_maternal_file <- "data/fixed_maternal_genotypes.csv"
if (file.exists(fixed_maternal_file)) {
  fixed_maternal_df <- read.table(fixed_maternal_file, sep = ",", header = TRUE, check.names = FALSE)
  
  # Format it into a population list compatible with our pipeline
  n_fixed <- nrow(fixed_maternal_df)
  
  # Ensure locus order matches SNPs_tested
  snp_dosages <- fixed_maternal_df[, SNPs_tested, drop = FALSE]
  
  fixed_genotypes <- lapply(seq_len(n_fixed), function(i) as.integer(snp_dosages[i, ]))
  
  fixed_maternal_pop <<- list(
    indiv_ID    = paste0("FIXED_", fixed_maternal_df$indiv_ID),
    genotypes   = fixed_genotypes,
    allele_freqs = original_init_population$allele_freqs,
    maternal_ID = rep(NA_integer_, n_fixed),
    paternal_ID = rep(NA_integer_, n_fixed),
    generation  = "FIXED",
    GEBV       = rep(0, n_fixed),
    error       = rep(0, n_fixed)
  )
  
  # GEBV for fixed maternal individuals is now calculated during population initialization.
  # We just retrieve it from the original_init_population object where it was appended.
  fixed_GEBVs <- original_init_population$GEBV[match(fixed_maternal_pop$indiv_ID, original_init_population$indiv_ID)]
  fixed_maternal_pop$GEBV = fixed_GEBVs
  
  cat(sprintf("Loaded %d fixed maternal individuals.\n", n_fixed))
} else {
  fixed_maternal_pop = NULL
}

## Helper: Subset a population to the top fraction by GEBV (can optionally lock in fixed maternal lines)
get_top_fraction <- function(pop, fraction, force_fixed_maternal = FALSE) {
  
  n_keep <- max(1, floor(length(pop$indiv_ID) * fraction))
  
  # Base selection: sort all individuals by GEBV (lower is better)
  idx <- order(pop$GEBV)[1:n_keep]
  
  selected_pop <- list(
    indiv_ID    = pop$indiv_ID[idx],
    genotypes   = pop$genotypes[idx],
    GEBV        = pop$GEBV[idx],
    allele_freqs = pop$allele_freqs,
    maternal_ID = pop$maternal_ID[idx],
    paternal_ID = pop$paternal_ID[idx],
    generation  = pop$generation,
    error       = pop$error[idx]
  )
  
  # If requested and available, replace the bottom of the selected pool with fixed individuals
  if (force_fixed_maternal && !is.null(fixed_maternal_pop)) {
    n_fixed <- length(fixed_maternal_pop$indiv_ID)
    
    # We can only displace up to n_keep individuals.
    n_replace <- min(n_fixed, n_keep)
    
    if (n_replace > 0) {
      cat(sprintf("  Displacing bottom %d individuals from the selected maternal pool with fixed lines.\n", n_replace))
      
      # The selected_pop is already sorted best-to-worst (indices 1 to n_keep)
      # We keep the best (1 : n_keep - n_replace), and append the fixed lines
      keep_idx <- seq_len(n_keep - n_replace)
      
      selected_pop$indiv_ID  <- c(selected_pop$indiv_ID[keep_idx], fixed_maternal_pop$indiv_ID[1:n_replace])
      selected_pop$genotypes <- c(selected_pop$genotypes[keep_idx], fixed_maternal_pop$genotypes[1:n_replace])
      selected_pop$GEBV      <- c(selected_pop$GEBV[keep_idx], fixed_maternal_pop$GEBV[1:n_replace])
      
      selected_pop$maternal_ID <- c(selected_pop$maternal_ID[keep_idx], fixed_maternal_pop$maternal_ID[1:n_replace])
      selected_pop$paternal_ID <- c(selected_pop$paternal_ID[keep_idx], fixed_maternal_pop$paternal_ID[1:n_replace])
      selected_pop$error       <- c(selected_pop$error[keep_idx], fixed_maternal_pop$error[1:n_replace])
    }
  }
  
  return(selected_pop)
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