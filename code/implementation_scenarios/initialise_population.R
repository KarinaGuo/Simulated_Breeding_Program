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
df_base <- read.csv(gt_datafile, sep = "\t", header = TRUE) 

SNP_AF_Histset <- read.csv("data/SNP_AF.csv", header = TRUE)

SNPs_tested <- intersect(intersect(SNPs_tested, SNP_AF_Histset$locus), df_base$rs.)

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
init_AF <- replicate(
  n        = population_size,
  expr     = rbinom(n    = length(allele_freqs),
                    size = 2,
                    prob = allele_freqs),
  simplify = FALSE
)

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
original_init_population <- list(
  indiv_ID    = indiv_ID,
  genotypes   = init_AF,        # list of N dosage vectors (length = n_SNPs)
  GEBV        = init_GEBV,      # numeric vector [N], rescaled to [0.001, 1]
  allele_freqs = allele_freqs,  # numeric vector [n_SNPs] — ALT freq from SNP_AF.csv
  maternal_ID = rep(NA_integer_, population_size),
  paternal_ID = rep(NA_integer_, population_size),
  generation  = "F0",
  error       = init_err
)

## ── Save ──────────────────────────────────────────────────────────────────────
saveRDS(original_init_population, file = "data/original_init_population.rds")
save.image(file = "data/original_data.rds")


cat("Initialisation complete.\n")
cat("  Population size :", population_size, "\n")
cat("  SNPs            :", length(SNPs_tested), "\n")
cat("  GEBV range      :", round(range(init_GEBV), 4), "\n")
cat("  AF range        :", round(range(allele_freqs), 4), "\n")
