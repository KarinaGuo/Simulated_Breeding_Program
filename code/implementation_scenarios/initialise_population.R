set.seed(12345)
rm(list = ls())
options(repos = c(CRAN = "https://cran.csiro.au/"))

library(tidyverse)
library(scales)

## ── Load functions ─────────────────────────────────────────────────────────────
## Note: all sources are relative to the code/ directory; set wd before running
source("code/gapit_functions_080425.txt")
source("code/GP_functs.R")
source("code/genotype_phenotype_v1.R")

## ── Parameters ────────────────────────────────────────────────────────────────
population_size <- 30

log_file            <- "data/pred_GAPIT_log.txt"
out_dir             <- "data"
gt_datafile         <- "data/Report-DMela25-10229/Report_DMela25-10229_RegularGenotyping/Report_DMela25-10229_GenotypingSamples_trainingconcat_sort.hapmap.hmp.txt"
input_phenodatafile <- "data/mq_phenotypes.csv"

## ── SNP list ──────────────────────────────────────────────────────────────────
SNP_list    <- read.csv("data/SNP_list.csv", header = TRUE)
SNPs_tested <- SNP_list[, 1]

indiv_ID <- seq_len(population_size)

## ── SNP metadata from training hapmap ─────────────────────────────────────────
df_base <- read.csv(gt_datafile, sep = "	", header = TRUE) 

SNP_AF_Histset <- read.csv("data/SNP_AF.csv", header = TRUE)

SNPs_tested <- intersect(intersect(SNPs_tested, SNP_AF_Histset$locus), df_base$rs.)

fixed_maternal_file <- "data/fixed_maternal_genotypes.csv"
if (file.exists(fixed_maternal_file)) {
  fixed_maternal_df <- read.table(fixed_maternal_file, sep = ",", header = TRUE, check.names = FALSE)
  SNPs_tested <- intersect(SNPs_tested, names(fixed_maternal_df))
}

uploaded_f0_file <- "data/uploaded_F0_genotypes.csv"
if (file.exists(uploaded_f0_file)) {
  uploaded_f0_df <- read.csv(uploaded_f0_file, check.names = FALSE)
  SNPs_tested <- intersect(SNPs_tested, names(uploaded_f0_df))
} else {
  uploaded_f0_df <- NULL
}


df_base <- df_base %>%
  filter(rs. %in% SNPs_tested)

SNP_AF_Histset <- SNP_AF_Histset %>% 
  filter(locus %in% SNPs_tested)


#----

hapmap_meta <- df_base[, 1:11]

training_hapmap <- df_base[, c(1, grep("S_", colnames(df_base)))]

snp_refalt_identity <- df_base %>%
  separate(
    col    = alleles,
    into   = c("REF", "ALT"),
    sep    = "/",
    fill   = "right",
    remove = FALSE
  ) %>%
  dplyr::select(rs., REF, ALT) %>% 
  filter(rs. %in% SNPs_tested)

## Reorder to match SNPs_tested order
SNPs_tested_ord           <- match(SNPs_tested, snp_refalt_identity$rs.)
stopifnot(all(!is.na(SNPs_tested_ord)))
snp_refalt_identity_reord <- snp_refalt_identity[SNPs_tested_ord, ]
ref_alleles               <- snp_refalt_identity_reord$REF
alt_alleles               <- snp_refalt_identity_reord$ALT

iupac_dict <- c(
  "AG" = "R", "GA" = "R",
  "CT" = "Y", "TC" = "Y",
  "GC" = "S", "CG" = "S",
  "AT" = "W", "TA" = "W",
  "GT" = "K", "TG" = "K",
  "AC" = "M", "CA" = "M"
)
het_codes <- iupac_dict[paste0(ref_alleles, alt_alleles)]

## ── Reference allele frequencies ──────────────────────────────────────────────
## SNP_AF.csv contains per-locus ALT allele frequencies estimated from the
## training population.  These are the REFERENCE allele frequencies for the
## whole project: they anchor kinship (VanRaden GRM), allelic diversity (He),
## and allelic capture comparisons across all generations.
##
## They are stored in `population$allele_freqs` so every downstream script
## can use the same baseline without re-loading the CSV.

## Ensure locus order matches SNPs_tested
af_ord       <- match(SNPs_tested, SNP_AF_Histset$locus)
stopifnot(all(!is.na(af_ord)))
allele_freqs <- SNP_AF_Histset$frequency[af_ord]   # ALT allele freq, one per SNP

## ── Simulate initial genotypes from allele frequencies ────────────────────────
## Each individual's genotype at locus l is drawn as Binomial(2, p_l),
## where p_l is the ALT allele frequency from SNP_AF.csv.
## This links each individual's genotype to the population-level AF baseline.

if (!is.null(uploaded_f0_df)) {
  population_size <- nrow(uploaded_f0_df)
  indiv_ID <- as.character(uploaded_f0_df$indiv_ID)
  
  # Ensure locus order matches the newly intersected SNPs_tested
  snp_dosages_f0 <- uploaded_f0_df[, SNPs_tested, drop = FALSE]
  init_AF <- lapply(seq_len(population_size), function(i) as.integer(snp_dosages_f0[i, ]))
  
  cat(sprintf("Loaded %d F0 individuals dynamically from uploaded file.\n", population_size))
} else {
  # Fallback to standard simulation if no custom F0 population is provided
  population_size <- 30
  indiv_ID <- as.character(seq_len(population_size))
  init_AF <- replicate(
    n        = population_size,
    expr     = rbinom(n    = length(allele_freqs),
                      size = 2,
                      prob = allele_freqs),
    simplify = FALSE
  )
  cat(sprintf("Simulated %d F0 individuals randomly.\n", population_size))
}
## ── Compute GEBV via GAPIT ────────────────────────────────────────────────────
## Phenotype_from_genotype_GAPIT is defined in genotype_phenotype_v1.R.
## It uses hapmap_meta, training_hapmap, ref_alleles, alt_alleles, het_codes,
## input_phenodatafile and out_dir from the enclosing environment.
init_GEBV <- Phenotype_from_genotype_GAPIT(
  individuals_GT = init_AF,
  SNPs_tested    = SNPs_tested
)

## Clip to [0.001, 1] as required by project rules
init_GEBV <- pmax(pmin(init_GEBV, 1), 0.001)

init_err <- rep(0, population_size)

## ── Build original_init_population ────────────────────────────────────────────
## allele_freqs is stored at the population level so that:
##   - summary_statistics.R uses it as the baseline for He and allelic capture
##   - calculate_kinship.R uses it to centre the GRM (VanRaden method)
##   - Phenotype_from_genotype_GAPIT uses the same locus-level AF when
##     constructing simulated hapmaps for subsequent generations

## Load Fixed Maternal Individuals if provided (optional)
## Expected format: A dataframe/CSV with columns indiv_ID, and SNP dosages (0,1,2) for all SNPs_tested
fixed_maternal_file <- "data/fixed_maternal_genotypes.csv"
if (file.exists(fixed_maternal_file)) {
  n_fixed <- nrow(fixed_maternal_df)
  
  # Ensure locus order matches SNPs_tested
  snp_dosages <- fixed_maternal_df[, SNPs_tested, drop = FALSE]
  fixed_genotypes <- lapply(seq_len(n_fixed), function(i) as.integer(snp_dosages[i, ]))
  
  fixed_indiv_ID <- paste0("FIXED_", fixed_maternal_df$indiv_ID)
  
  # Calculate GEBV for fixed maternal individuals
  fixed_GEBV <- pmax(pmin(Phenotype_from_genotype_GAPIT(fixed_genotypes, SNPs_tested), 1), 0.001)
  
  # Combine with random initial population
  combined_indiv_ID <- as.character(c(indiv_ID, fixed_indiv_ID))
  combined_AF <- c(init_AF, fixed_genotypes)
  combined_GEBV <- c(init_GEBV, fixed_GEBV)
  combined_population_size <- population_size + n_fixed
  
  original_init_population <- list(
    indiv_ID    = combined_indiv_ID,
    genotypes   = combined_AF,        # list of N dosage vectors (length = n_SNPs)
    GEBV        = combined_GEBV,      # numeric vector [N], rescaled to [0.001, 1]
    allele_freqs = allele_freqs,  # numeric vector [n_SNPs] — ALT freq from SNP_AF.csv
    maternal_ID = rep(NA_character_, combined_population_size),
    paternal_ID = rep(NA_character_, combined_population_size),
    generation  = "F0",
    error       = rep(0, combined_population_size)
  )
  cat(sprintf("Loaded %d fixed maternal individuals into the initial population.\n", n_fixed))
} else {
  original_init_population <- list(
    indiv_ID    = as.character(indiv_ID),
    genotypes   = init_AF,        # list of N dosage vectors (length = n_SNPs)
    GEBV        = init_GEBV,      # numeric vector [N], rescaled to [0.001, 1]
    allele_freqs = allele_freqs,  # numeric vector [n_SNPs] — ALT freq from SNP_AF.csv
    maternal_ID = rep(NA_character_, population_size),
    paternal_ID = rep(NA_character_, population_size),
    generation  = "F0",
    error       = init_err
  )
}


## ── Save ──────────────────────────────────────────────────────────────────────
saveRDS(original_init_population, file = "data/original_init_population.rds")
save.image(file = "data/original_data.rds")


cat("Initialisation complete.\n")
cat("  Population size :", population_size, "\n")
cat("  SNPs            :", length(SNPs_tested), "\n")
cat("  GEBV range      :", round(range(init_GEBV), 4), "\n")
cat("  AF range        :", round(range(allele_freqs), 4), "\n")