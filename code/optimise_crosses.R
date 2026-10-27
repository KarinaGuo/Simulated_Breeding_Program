## optimise_crosses.R
##
## Kinship + GEBV minimiser for breeding cross selection.
## Adapted from Jason Bragg's OptGenMix package
## (OptGenMix-master/OptGenMix-master/R/).
##
## Purpose (Scenarios 3, 4, 7):
##   Given a pool of candidate parents, select the PATERNAL parent for each
##   maternal individual by minimising a combined score:
##
##     score(mat, pat) = (1 - gebv_weight) * kinship(mat, pat)
##                     +      gebv_weight  * GEBV(pat)
##
##   Lower score = better (project rule: lower GEBV is optimal).
##
##   Two strategies are provided:
##     1. greedy_select_paternal()  — fast, deterministic, greedy (good default)
##     2. sa_select_paternal()      — simulated annealing over all cross
##                                    assignments simultaneously (global search,
##                                    mirrors OptGenMix's optimize_single_objective)
##
## Key differences from OptGenMix:
##   * Input is the project population list (not a raw genotype / similarity matrix)
##   * Objective is minimisation (lower GEBV = better), not maximisation
##   * Assignment is one paternal per maternal (1-to-1 with repeats), not weights
##   * No genotype diversity measures (Nei / Shannon) — GEBV + kinship only
##   * generate_offspring() in optimise_grid_placement.R handles Mendelian crosses
##
## Usage:
##   source("calculate_kinship.R")
##   source("optimise_crosses.R")
##
##   K <- calc_kinship(pop, allele_freqs = original_init_population$allele_freqs)
##
##   ## Greedy (fast, good for most scenarios):
##   plan <- greedy_select_paternal(maternal_pop, paternal_pool, K,
##                                  gebv_weight = 0.5)
##
##   ## SA (global, slower — use for Scenario 7 with all individuals):
##   plan <- sa_select_paternal(maternal_pop, paternal_pool, K,
##                              gebv_weight = 0.5, max_steps = 20000)
##
##   ## Then generate offspring (defined in optimise_grid_placement.R):
##   offspring <- generate_offspring(plan, parent_pop, generation = "F1")

## ═══════════════════════════════════════════════════════════════════════════════
## PAIR SCORING
## ═══════════════════════════════════════════════════════════════════════════════

## ── Score a single cross ──────────────────────────────────────────────────────
## @param mat_id      Character. Maternal individual ID.
## @param pat_id      Character. Paternal individual ID.
## @param kinship_mat N×N GRM matrix (named by indiv_ID)
## @param gebv_vec    Named numeric vector of GEBV scores (names = indiv_ID)
## @param gebv_weight Numeric in [0,1].  0 = kinship only; 1 = GEBV only.
## @return Single numeric score (lower = better)

score_cross <- function(mat_id, pat_id, kinship_mat, gebv_vec,
                         gebv_weight = 0.5) {
  kin  <- kinship_mat[as.character(mat_id), as.character(pat_id)]
  gebv <- gebv_vec[as.character(pat_id)]
  (1 - gebv_weight) * kin + gebv_weight * gebv
}

## ── Score an entire crossing plan ─────────────────────────────────────────────
## @param plan        data.frame with columns maternal_ID, paternal_ID
## @param kinship_mat N×N GRM matrix
## @param gebv_vec    Named GEBV vector
## @param gebv_weight Numeric weight
## @return Mean cross score across all maternal–paternal pairs (lower = better)

score_plan <- function(plan, kinship_mat, gebv_vec, gebv_weight = 0.5) {
  scores <- mapply(
    score_cross,
    mat_id      = as.character(plan$maternal_ID),
    pat_id      = as.character(plan$paternal_ID),
    MoreArgs    = list(kinship_mat = kinship_mat,
                       gebv_vec    = gebv_vec,
                       gebv_weight = gebv_weight)
  )
  mean(scores, na.rm = TRUE)
}

## ── Full cross-score matrix ───────────────────────────────────────────────────
## Builds an n_maternal × n_paternal matrix of pair scores.
## Used by both greedy and SA methods.
##
## @param maternal_ids Character vector of maternal IDs
## @param paternal_ids Character vector of paternal IDs (candidate pool)
## @param kinship_mat  N×N GRM
## @param gebv_vec     Named GEBV vector
## @param gebv_weight  Numeric weight
## @return Named numeric matrix (maternal × paternal), lower = better

build_score_matrix <- function(maternal_ids, paternal_ids,
                                kinship_mat, gebv_vec,
                                gebv_weight = 0.5) {
  mat_ids <- as.character(maternal_ids)
  pat_ids <- as.character(paternal_ids)

  ## Kinship sub-matrix (maternal rows × paternal cols)
  kin_sub  <- kinship_mat[mat_ids, pat_ids, drop = FALSE]

  ## GEBV row (repeated for each maternal)
  gebv_row <- matrix(gebv_vec[pat_ids],
                     nrow = length(mat_ids),
                     ncol = length(pat_ids),
                     byrow = TRUE)

  score_mat <- (1 - gebv_weight) * kin_sub + gebv_weight * gebv_row
  rownames(score_mat) <- mat_ids
  colnames(score_mat) <- pat_ids
  score_mat
}

## ═══════════════════════════════════════════════════════════════════════════════
## STRATEGY 1 — GREEDY SELECTION
## ═══════════════════════════════════════════════════════════════════════════════

## ── Greedy paternal selection ─────────────────────────────────────────────────
## For each maternal individual (in order of worst GEBV first, so the hardest
## individuals get the best available paternal mate), assigns the paternal with
## the lowest combined score.  Repeats are allowed in the paternal pool.
##
## This is fast and performs well when the paternal pool is larger than the
## maternal set (Scenarios 3, 4: paternal drawn from same top-third/half pool).
##
## @param maternal_pop  Population list for maternal individuals (sorted by GEBV)
## @param paternal_pool Population list for paternal candidates
## @param kinship_mat   Full N×N GRM (must include both maternal and paternal IDs)
## @param gebv_weight   Weight on GEBV component (default 0.5)
## @param target_n      Target total offspring (~100); distributes equally across
##                      maternal parents (project rule: equal offspring per individual)
## @param allow_selfing Logical. If FALSE (default), maternal == paternal is forbidden.
##
## @return data.frame with columns: maternal_ID, paternal_ID, n_offspring, cross_score

greedy_select_paternal <- function(maternal_pop, paternal_pool,
                                   kinship_mat,
                                   gebv_weight   = 0.5,
                                   target_n      = 100,
                                   allow_selfing = FALSE) {

  mat_ids <- as.character(maternal_pop$indiv_ID)
  pat_ids <- as.character(paternal_pool$indiv_ID)

  gebv_vec <- c(
    setNames(maternal_pop$GEBV, mat_ids),
    setNames(paternal_pool$GEBV, pat_ids)
  )
  gebv_vec <- gebv_vec[!duplicated(names(gebv_vec))]

  ## Build pair score matrix
  score_mat <- build_score_matrix(mat_ids, pat_ids,
                                  kinship_mat, gebv_vec, gebv_weight)

  n_mat <- length(mat_ids)

  ## Offspring distribution
  n_each  <- floor(target_n / n_mat)
  n_extra <- target_n - n_each * n_mat

  ## Process maternals in order of worst (highest) GEBV first — they benefit
  ## most from a carefully chosen paternal partner
  order_mat <- order(gebv_vec[mat_ids], decreasing = TRUE)

  plan <- vector("list", n_mat)

  for (rank in seq_len(n_mat)) {
    i      <- order_mat[rank]
    mat_id <- mat_ids[i]

    ## Scores for this maternal against all paternals
    pair_scores <- score_mat[mat_id, ]

    ## Forbid selfing if requested
    if (!allow_selfing && mat_id %in% pat_ids) {
      pair_scores[mat_id] <- Inf
    }

    ## Pick the paternal with the lowest combined score
    best_pat <- names(which.min(pair_scores))
    best_score <- pair_scores[best_pat]

    n_off <- n_each + if (rank <= n_extra) 1L else 0L

    plan[[i]] <- data.frame(
      maternal_ID = mat_id,
      paternal_ID = best_pat,
      n_offspring = n_off,
      cross_score = best_score,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, plan)
}

## ═══════════════════════════════════════════════════════════════════════════════
## STRATEGY 2 — SIMULATED ANNEALING (GLOBAL SEARCH)
## ═══════════════════════════════════════════════════════════════════════════════

## ── SA paternal selection ─────────────────────────────────────────────────────
## Optimises the full assignment of paternal partners simultaneously via
## simulated annealing, inspired by OptGenMix's optimize_single_objective.
##
## At each SA step, one maternal individual's paternal partner is changed to a
## randomly drawn alternative from the paternal pool.  The move is accepted if
## it lowers the mean plan score, or with probability exp(-delta/T) otherwise.
##
## Use this when:
##   - You want a global optimum rather than a greedy solution
##   - The paternal pool is the same as the maternal set (Scenario 7: all individuals)
##     because greedy can get trapped by early high-kinship assignments
##
## @param maternal_pop   Population list for maternal individuals
## @param paternal_pool  Population list for paternal candidates (can = maternal_pop)
## @param kinship_mat    Full N×N GRM
## @param gebv_weight    Weight on GEBV (default 0.5)
## @param target_n       Target offspring total
## @param allow_selfing  Logical. Default FALSE.
## @param max_steps      SA chain length (default 20000)
## @param max_t          Starting temperature (default 0.01)
## @param seed           Optional integer seed
##
## @return data.frame with columns: maternal_ID, paternal_ID, n_offspring,
##         cross_score, plus attribute "chain" (step, score_current, accepted)

sa_select_paternal <- function(maternal_pop, paternal_pool,
                                kinship_mat,
                                gebv_weight   = 0.5,
                                target_n      = 100,
                                allow_selfing = FALSE,
                                max_steps     = 20000,
                                max_t         = 0.01,
                                seed          = NULL) {

  if (!is.null(seed)) set.seed(seed)

  mat_ids <- as.character(maternal_pop$indiv_ID)
  pat_ids <- as.character(paternal_pool$indiv_ID)
  n_mat   <- length(mat_ids)

  gebv_vec <- c(
    setNames(maternal_pop$GEBV, mat_ids),
    setNames(paternal_pool$GEBV, pat_ids)
  )
  gebv_vec <- gebv_vec[!duplicated(names(gebv_vec))]

  ## Pre-build score matrix for fast lookup
  score_mat <- build_score_matrix(mat_ids, pat_ids,
                                  kinship_mat, gebv_vec, gebv_weight)

  ## Offspring distribution
  n_each  <- floor(target_n / n_mat)
  n_extra <- target_n - n_each * n_mat

  ## ── Initialise: random assignment ────────────────────────────────────────────
  ## For each maternal, sample a paternal at random (no selfing if requested)
  current_pat <- character(n_mat)
  names(current_pat) <- mat_ids

  for (i in seq_len(n_mat)) {
    candidates <- pat_ids
    if (!allow_selfing) candidates <- setdiff(candidates, mat_ids[i])
    current_pat[i] <- sample(candidates, 1L)
  }

  ## Current mean plan score (look up from score matrix)
  cur_score <- mean(score_mat[cbind(mat_ids, current_pat)], na.rm = TRUE)

  ## Chain storage (compact: just score and acceptance)
  chain <- data.frame(
    step          = seq_len(max_steps),
    score_current = NA_real_,
    accepted      = FALSE
  )

  cat(sprintf(
    "SA cross selection: %d maternal × %d paternal, %d steps, gebv_weight = %.2f\n",
    n_mat, length(pat_ids), max_steps, gebv_weight
  ))
  cat(sprintf("  Initial mean score: %.6f\n", cur_score))

  ## ── SA main loop ─────────────────────────────────────────────────────────────
  for (s in seq_len(max_steps)) {

    ## Linear cooling (same schedule as OptGenMix / Jason Bragg scripts)
    t <- (1 - (s - 1) / max_steps) * max_t

    ## Propose: change one randomly chosen maternal's paternal partner
    i_change <- sample.int(n_mat, 1L)
    mat_id   <- mat_ids[i_change]

    candidates <- pat_ids
    if (!allow_selfing) candidates <- setdiff(candidates, mat_id)

    ## Exclude current paternal to ensure a real change
    candidates <- setdiff(candidates, current_pat[i_change])
    if (length(candidates) == 0L) {
      chain$score_current[s] <- cur_score
      next
    }

    new_pat <- sample(candidates, 1L)

    ## Score delta for the changed pair only (efficient: don't recompute whole plan)
    old_pair_score <- score_mat[mat_id, current_pat[i_change]]
    new_pair_score <- score_mat[mat_id, new_pat]

    ## Change in mean plan score
    delta_mean <- (new_pair_score - old_pair_score) / n_mat
    new_score  <- cur_score + delta_mean

    ## Accept / reject
    accepted <- FALSE
    if (new_score < cur_score) {
      accepted <- TRUE
    } else if (t > 0) {
      if (runif(1) < exp(-delta_mean / t)) accepted <- TRUE
    }

    if (accepted) {
      current_pat[i_change] <- new_pat
      cur_score             <- new_score
    }

    chain$score_current[s] <- cur_score
    chain$accepted[s]      <- accepted

    if (s %% max(1L, floor(max_steps / 10L)) == 0L) {
      cat(sprintf("  Step %6d / %d — mean score: %.6f\n", s, max_steps, cur_score))
    }
  }

  cat(sprintf("  Final mean score : %.6f\n", cur_score))

  ## ── Assemble output plan ──────────────────────────────────────────────────────
  plan <- data.frame(
    maternal_ID = mat_ids,
    paternal_ID = as.character(current_pat),
    n_offspring = n_each + c(rep(1L, n_extra), rep(0L, n_mat - n_extra)),
    cross_score = score_mat[cbind(mat_ids, as.character(current_pat))],
    stringsAsFactors = FALSE
  )
  attr(plan, "chain") <- chain
  plan
}

## ═══════════════════════════════════════════════════════════════════════════════
## DIAGNOSTICS
## ═══════════════════════════════════════════════════════════════════════════════

## ── Summarise a crossing plan ─────────────────────────────────────────────────
## Prints key statistics about a plan returned by greedy or SA selection.
##
## @param plan        data.frame from greedy_select_paternal() or sa_select_paternal()
## @param kinship_mat N×N GRM
## @param gebv_vec    Named GEBV vector

summarise_plan <- function(plan, kinship_mat, gebv_vec) {
  kin_vals  <- mapply(function(m, p) kinship_mat[as.character(m), as.character(p)],
                      plan$maternal_ID, plan$paternal_ID)
  gebv_vals <- gebv_vec[as.character(plan$paternal_ID)]

  cat("── Crossing plan summary ──────────────────────────────\n")
  cat(sprintf("  Pairs (maternal × paternal) : %d\n",       nrow(plan)))
  cat(sprintf("  Total offspring target      : %d\n",       sum(plan$n_offspring)))
  cat(sprintf("  Mean pairwise kinship       : %.4f\n",     mean(kin_vals, na.rm = TRUE)))
  cat(sprintf("  Mean paternal GEBV          : %.4f\n",     mean(gebv_vals, na.rm = TRUE)))
  cat(sprintf("  Mean cross score            : %.4f\n",     mean(plan$cross_score, na.rm = TRUE)))
  cat(sprintf("  Unique paternal parents     : %d / %d\n",
              length(unique(plan$paternal_ID)), nrow(plan)))
  cat("───────────────────────────────────────────────────────\n")
}

## ── Plot SA convergence chain ─────────────────────────────────────────────────
## @param plan  Output from sa_select_paternal() (must have "chain" attribute)

plot_cross_chain <- function(plan) {
  chain <- attr(plan, "chain")
  if (is.null(chain)) {
    message("plot_cross_chain: no chain attribute (only sa_select_paternal produces one).")
    return(invisible(NULL))
  }
  plot(
    chain$step, chain$score_current,
    type = "l", col = "steelblue", lwd = 1.5,
    xlab = "SA step",
    ylab = "Mean cross score (lower = better)",
    main = "Cross selection SA convergence"
  )
}
