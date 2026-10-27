## calculate_kinship.R
##
## Calculates a genomic relationship (kinship) matrix from population genotypes.
## Uses the VanRaden (2008) method — the same approach underlying gBLUP in GAPIT.
##
## Dependencies: none (base R only)
##
## Usage:
##   source("calculate_kinship.R")
##
##   ## Using the reference allele frequencies from the initial population:
##   K <- calc_kinship(current_pop, allele_freqs = original_init_population$allele_freqs)
##
## The VanRaden GRM centres individuals relative to the F0 reference allele
## frequencies stored in original_init_population$allele_freqs.  This is
## essential for comparability across generations:  an F2 individual with the
## same genotype as an F0 individual will have the same row/column in K
## regardless of the generation's observed allele frequencies.
##
## Reference:
##   VanRaden, P.M. (2008). Efficient methods to compute genomic predictions.
##   Journal of Dairy Science, 91(11), 4414–4423.

## ── VanRaden GRM ──────────────────────────────────────────────────────────────
##
## Z = M - 2P   where M[i,l] = dosage of individual i at locus l (0/1/2)
##              and P is a matrix where each column is 2 * p_l
##              with p_l the REFERENCE (F0) ALT allele frequency from SNP_AF.csv
##
## G = Z Z' / (2 * sum_l p_l * (1 - p_l))
##
## @param pop          Population list with $genotypes (list of dosage vectors)
## @param allele_freqs Numeric vector of reference ALT allele frequencies,
##                     one per SNP.  Should be original_init_population$allele_freqs.
##                     If NULL, frequencies are estimated from pop$genotypes
##                     (acceptable for exploratory use but not for cross-generation
##                     comparison).
##
## @return Symmetric N x N genomic relationship matrix (named by indiv_ID)

calc_kinship <- function(pop, allele_freqs = NULL) {

  stopifnot(
    is.list(pop),
    "genotypes" %in% names(pop)
  )

  ## Build N x L dosage matrix  (rows = individuals, cols = SNPs)
  M <- do.call(rbind, lapply(pop$genotypes, as.numeric))   # N x L

  n_indiv <- nrow(M)
  n_loci  <- ncol(M)

  ## Allele frequencies
  if (is.null(allele_freqs)) {
    ## Estimate from current population dosages.
    ## NOTE: for cross-generation comparisons use original_init_population$allele_freqs
    warning("calc_kinship: allele_freqs not supplied — estimating from current population. ",
            "Pass original_init_population$allele_freqs for consistent cross-generation GRM.")
    p <- colMeans(M, na.rm = TRUE) / 2
  } else {
    stopifnot(length(allele_freqs) == n_loci)
    p <- allele_freqs
  }

  ## Centering matrix: subtract 2*p from each column
  P <- matrix(2 * p, nrow = n_indiv, ncol = n_loci, byrow = TRUE)
  Z <- M - P

  ## Scaling factor
  scale_factor <- 2 * sum(p * (1 - p), na.rm = TRUE)

  if (scale_factor == 0) {
    stop("calc_kinship: scaling factor is 0 — all SNPs may be monomorphic.")
  }

  ## GRM
  G <- (Z %*% t(Z)) / scale_factor

  ## Name rows/cols by individual ID if available
  if (!is.null(pop$indiv_ID)) {
    rownames(G) <- colnames(G) <- as.character(pop$indiv_ID)
  }

  G
}

## ── Subset kinship matrix ─────────────────────────────────────────────────────
## Return the kinship sub-matrix for a subset of individual IDs.
##
## @param K       Full N x N kinship matrix (named rows/cols)
## @param id_vec  Character or integer vector of individual IDs to keep

subset_kinship <- function(K, id_vec) {
  id_vec <- as.character(id_vec)
  K[id_vec, id_vec, drop = FALSE]
}

## ── Mean pairwise kinship ─────────────────────────────────────────────────────
## Returns the mean off-diagonal kinship (average relatedness in population).

mean_pairwise_kinship <- function(K) {
  off <- K[lower.tri(K)]
  mean(off, na.rm = TRUE)
}

## ── Kinship between two individuals ──────────────────────────────────────────
## Convenience wrapper used in pairing optimisation.
##
## @param K    Kinship matrix
## @param id_i ID of individual i (must match rownames of K)
## @param id_j ID of individual j

pairwise_kinship <- function(K, id_i, id_j) {
  K[as.character(id_i), as.character(id_j)]
}
