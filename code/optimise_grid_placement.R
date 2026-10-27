## optimise_grid_placement.R
##
## Grid-based spatial optimisation for the breeding project.
## Adapted from Jason Bragg's spatial-genetic distance scripts
## (optimisation_scripts.R / optimisation_workflow.r).
##
## Purpose (Scenario 8):
##   1. Build the most compact (square-ish) grid that fits N individuals.
##   2. Place individuals in the grid to minimise the combined objective:
##        α * mean(kinship_ij / dist_ij) + (1-α) * mean(GEBV_ij / dist_ij)
##      where GEBV_ij = GEBV_i * GEBV_j (product — low when either parent is good).
##      This separates related individuals spatially AND biases frequent-crossers
##      toward individuals with low (desirable) GEBV.
##   3. Generate a crossing plan weighted by inverse spatial distance:
##      nearby individuals are more likely to cross.
##
## Key design choices vs the original scripts:
##   * Vectorised pair scoring — no slow nested loops in the SA chain
##   * Combined kinship + GEBV objective with a tuneable weight parameter
##   * Targeted swap (force worst-scoring individual) + random swap, alternating
##   * Output integrates with the project population list format
##
## Usage:
##   source("calculate_kinship.R")
##   source("optimise_grid_placement.R")
##
##   K      <- calc_kinship(pop, allele_freqs = original_init_population$allele_freqs)
##   result <- optimise_grid_placement(pop, K,
##                                     gebv_weight = 0.5,
##                                     max_steps   = 50000,
##                                     max_t       = 0.00001,
##                                     seed        = 42)
##   crosses <- grid_crosses_from_solution(result, pop, target_n = 100)

## ═══════════════════════════════════════════════════════════════════════════════
## GRID CONSTRUCTION
## ═══════════════════════════════════════════════════════════════════════════════

## ── Build square(ish) grid ────────────────────────────────────────────────────
## Finds the most compact rectangular grid (minimising perimeter = 2*(x+y))
## that fits exactly n_indiv cells.  For any N:
##   x = ceiling(sqrt(N)),  y = ceiling(N / x)
## gives x*y >= N with the smallest perimeter among integer rectangles.
##
## If x*y > N the last (x*y - N) cells (bottom-right) are dropped, keeping the
## grid contiguous.
##
## @param n_indiv  Integer. Number of individuals / grid sites needed.
## @return data.frame with columns: location_name, x_pos, y_pos

build_grid <- function(n_indiv) {
  x <- ceiling(sqrt(n_indiv))
  y <- ceiling(n_indiv / x)
  
  grid_coords <- expand.grid(x_pos = seq_len(x), y_pos = seq_len(y))
  
  ## Drop excess cells from the end
  if (nrow(grid_coords) > n_indiv) {
    grid_coords <- grid_coords[seq_len(n_indiv), ]
  }
  
  grid_coords$location_name <- paste0("L", seq_len(nrow(grid_coords)))
  rownames(grid_coords)     <- grid_coords$location_name
  grid_coords[, c("location_name", "x_pos", "y_pos")]
}

## ── Spatial distance matrix ───────────────────────────────────────────────────
## Euclidean distances between all pairs of grid positions.
## @param grid  data.frame from build_grid()
## @return Named numeric matrix (n_sites × n_sites)

build_spatial_distance <- function(grid) {
  sd <- as.matrix(dist(grid[, c("x_pos", "y_pos")]))
  rownames(sd) <- colnames(sd) <- grid$location_name
  sd
}

## ── Grid summary: perimeter and density ───────────────────────────────────────
## Prints a short description of the grid dimensions and compactness.
## @param grid  data.frame from build_grid()

describe_grid <- function(grid) {
  x_max <- max(grid$x_pos)
  y_max <- max(grid$y_pos)
  n     <- nrow(grid)
  cat(sprintf(
    "Grid: %d × %d = %d cells (%d individuals), perimeter = %d\n",
    x_max, y_max, x_max * y_max, n, 2L * (x_max + y_max)
  ))
}

## ═══════════════════════════════════════════════════════════════════════════════
## SCORING (VECTORISED)
## ═══════════════════════════════════════════════════════════════════════════════

## ── Precompute pair indices ───────────────────────────────────────────────────
## Returns a 2-column integer matrix of all upper-triangle pair indices (i < j).
## Computed once before the SA loop, not inside it.
## @param n  Number of locations
## @return   Integer matrix, nrow = n*(n-1)/2, ncol = 2

pair_indices <- function(n) {
  idx <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)
  idx  # col 1 = i, col 2 = j  (row indices)
}

## ── Vectorised pair score ─────────────────────────────────────────────────────
## For each pair (i,j) computes:
##   kinship_component = kinship[id_i, id_j] / spatial_dist[loc_i, loc_j]
##   gebv_component    = (GEBV[id_i] * GEBV[id_j]) / spatial_dist[loc_i, loc_j]
##
## Combined score:  α * mean(kinship_component) + (1-α) * mean(gebv_component)
## We MINIMISE this: small score means related/high-GEBV pairs are far apart.
##
## @param assignment   Named character vector (location_name → indiv_ID)
## @param kinship_mat  N×N matrix, named by indiv_ID
## @param gebv_vec     Named numeric vector of GEBV scores, names = indiv_ID
## @param spatial_dist n_sites×n_sites matrix, named by location_name
## @param pairs        Integer matrix from pair_indices()
## @param gebv_weight  Numeric in [0,1]. Weight on GEBV component (default 0.5)
## @return Single numeric score (lower = better)

score_assignment_vec <- function(assignment, kinship_mat, gebv_vec,
                                 spatial_dist, pairs, gebv_weight = 0.5) {
  loc_names <- names(assignment)
  ids       <- as.character(assignment)
  
  ## Spatial distances for all pairs (preindexed)
  sp <- spatial_dist[pairs]                        # n_pairs numeric vector
  ok <- sp > 0                                     # exclude same-cell (dist=0)
  
  if (!any(ok)) return(Inf)
  
  ## Kinship for all pairs
  kin_pairs <- kinship_mat[ids[pairs[, 1]], ids[pairs[, 2]]]
  if (is.matrix(kin_pairs)) kin_pairs <- diag(kin_pairs)  # safety
  
  ## GEBV product for all pairs  (lower GEBV is better; product is small when
  ## either or both parents have low GEBV — encourages central placement of
  ## good individuals who then cross frequently with neighbours)
  g_pairs <- gebv_vec[ids[pairs[, 1]]] * gebv_vec[ids[pairs[, 2]]]
  
  ## Weighted combined score
  kin_score  <- mean(kin_pairs[ok] / sp[ok],  na.rm = TRUE)
  gebv_score <- mean(g_pairs[ok]   / sp[ok],  na.rm = TRUE)
  
  (1 - gebv_weight) * kin_score + gebv_weight * gebv_score
}

## ═══════════════════════════════════════════════════════════════════════════════
## SIMULATED ANNEALING
## ═══════════════════════════════════════════════════════════════════════════════

## ── Propose a swap ────────────────────────────────────────────────────────────
## Swaps two individuals between grid locations.
## If `force_idx` is given, one of the two swap positions is fixed (used to
## target the worst-scoring individual — see original optimize_locations logic).
##
## @param assignment Named character vector (location → indiv_ID)
## @param force_idx  Optional integer index to include in the swap
## @return Updated assignment vector

propose_switch <- function(assignment, force_idx = NULL) {
  n <- length(assignment)
  if (is.null(force_idx)) {
    idx <- sample.int(n, 2L, replace = FALSE)
  } else {
    j   <- sample(setdiff(seq_len(n), force_idx), 1L)
    idx <- c(force_idx, j)
  }
  new_a        <- assignment
  new_a[idx[1]] <- assignment[idx[2]]
  new_a[idx[2]] <- assignment[idx[1]]
  new_a
}

## ── Find worst-scoring location ───────────────────────────────────────────────
## Returns the index of the location involved in the highest-scoring (worst) pair.
## This is used for targeted swaps (as in Jason Bragg's optimize_locations).
##
## @param assignment   Current assignment
## @param kinship_mat  Kinship matrix
## @param gebv_vec     Named GEBV vector
## @param spatial_dist Spatial distance matrix
## @param pairs        pair_indices() matrix
## @param gebv_weight  Weight on GEBV component
## @return Integer index into assignment

worst_location_idx <- function(assignment, kinship_mat, gebv_vec,
                               spatial_dist, pairs, gebv_weight) {
  ids      <- as.character(assignment)
  sp       <- spatial_dist[pairs]
  ok       <- sp > 0
  
  kin_pairs  <- kinship_mat[ids[pairs[, 1]], ids[pairs[, 2]]]
  if (is.matrix(kin_pairs)) kin_pairs <- diag(kin_pairs)
  
  g_pairs    <- gebv_vec[ids[pairs[, 1]]] * gebv_vec[ids[pairs[, 2]]]
  
  pairwise_score        <- rep(NA_real_, nrow(pairs))
  pairwise_score[ok]    <- (1 - gebv_weight) * kin_pairs[ok] / sp[ok] +
    gebv_weight       * g_pairs[ok]  / sp[ok]
  
  worst_pair <- which.max(pairwise_score)
  ## Randomly choose one of the two individuals in the worst pair
  sample(c(pairs[worst_pair, 1], pairs[worst_pair, 2]), 1L)
}

## ── Main SA optimiser ─────────────────────────────────────────────────────────
##
## @param pop          Population list with $indiv_ID and $GEBV
## @param kinship_mat  N×N GRM from calc_kinship(), named by indiv_ID
## @param gebv_weight  Weight on GEBV component of objective (default 0.5).
##                     0 = kinship only; 1 = GEBV only; 0.5 = equal weight.
## @param max_steps    SA chain length (default 50000)
## @param max_t        Starting temperature (default 0.00001)
## @param targeted_freq Fraction of steps that use a targeted (worst-pair) swap
##                       vs a random swap. Default 0.5.
## @param seed         Optional integer random seed
##
## @return Named list:
##   $solution     Named character vector: location_name → indiv_ID (final)
##   $grid         data.frame from build_grid()
##   $chain        data.frame (step, score_current, score_proposed, accepted)
##   $final_score  Numeric, objective value at solution
##   $gebv_weight  The weight used

optimise_grid_placement <- function(pop, kinship_mat,
                                    gebv_weight    = 0.5,
                                    max_steps      = 50000,
                                    max_t          = 0.00001,
                                    targeted_freq  = 0.5,
                                    seed           = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  stopifnot(
    is.list(pop),
    all(c("indiv_ID", "GEBV") %in% names(pop)),
    is.matrix(kinship_mat),
    nrow(kinship_mat) == length(pop$indiv_ID),
    gebv_weight >= 0, gebv_weight <= 1
  )
  
  n_indiv   <- length(pop$indiv_ID)
  indiv_ids <- as.character(pop$indiv_ID)
  
  ## Named GEBV vector for fast lookup
  gebv_vec <- setNames(pop$GEBV, indiv_ids)
  
  ## Ensure kinship matrix rows/cols match indiv_ids
  stopifnot(all(indiv_ids %in% rownames(kinship_mat)))
  
  ## ── Grid setup ──────────────────────────────────────────────────────────────
  grid         <- build_grid(n_indiv)
  describe_grid(grid)
  spatial_dist <- build_spatial_distance(grid)
  loc_names    <- grid$location_name
  n_locs       <- length(loc_names)
  
  ## Precompute pair index matrix (used for vectorised scoring throughout)
  pairs <- pair_indices(n_locs)
  
  ## ── Initial random assignment ────────────────────────────────────────────────
  assignment        <- sample(indiv_ids)
  names(assignment) <- loc_names
  
  cur_score <- score_assignment_vec(assignment, kinship_mat, gebv_vec,
                                    spatial_dist, pairs, gebv_weight)
  
  ## ── Chain storage ────────────────────────────────────────────────────────────
  chain <- data.frame(
    step           = seq_len(max_steps),
    score_current  = NA_real_,
    score_proposed = NA_real_,
    accepted       = FALSE
  )
  
  cat(sprintf(
    "SA grid optimisation: %d individuals, %d steps, gebv_weight = %.2f\n",
    n_indiv, max_steps, gebv_weight
  ))
  cat(sprintf("  Initial score: %.6f\n", cur_score))
  
  ## ── SA main loop ─────────────────────────────────────────────────────────────
  for (s in seq_len(max_steps)) {
    
    ## Linear cooling schedule (identical to Jason Bragg's optimize_locations_sa)
    t <- (1 - (s - 1) / max_steps) * max_t
    
    ## Alternate between targeted and random swaps
    if (runif(1) < targeted_freq) {
      widx          <- worst_location_idx(assignment, kinship_mat, gebv_vec,
                                          spatial_dist, pairs, gebv_weight)
      new_assignment <- propose_switch(assignment, force_idx = widx)
    } else {
      new_assignment <- propose_switch(assignment)
    }
    
    new_score <- score_assignment_vec(new_assignment, kinship_mat, gebv_vec,
                                      spatial_dist, pairs, gebv_weight)
    
    ## Accept / reject (Metropolis criterion)
    accepted <- FALSE
    if (new_score < cur_score) {
      accepted <- TRUE
    } else if (t > 0) {
      delta    <- new_score - cur_score
      if (runif(1) < exp(-delta / t)) accepted <- TRUE
    }
    
    if (accepted) {
      assignment <- new_assignment
      cur_score  <- new_score
    }
    
    chain$score_current[s]  <- cur_score
    chain$score_proposed[s] <- new_score
    chain$accepted[s]       <- accepted
    
    if (s %% max(1L, floor(max_steps / 10L)) == 0L) {
      cat(sprintf("  Step %6d / %d — score: %.6f\n", s, max_steps, cur_score))
    }
  }
  
  cat(sprintf("  Final score  : %.6f\n", cur_score))
  
  list(
    solution    = assignment,
    grid        = grid,
    chain       = chain,
    final_score = cur_score,
    gebv_weight = gebv_weight
  )
}

## ═══════════════════════════════════════════════════════════════════════════════
## CROSSING PLAN
## ═══════════════════════════════════════════════════════════════════════════════

## ── Generate crossing plan weighted by inverse spatial distance ───────────────
## For Scenario 8: each individual is the maternal parent exactly once.
## The paternal partner is sampled with probability ∝ 1/distance from the
## maternal individual's grid location (nearby individuals cross more often).
##
## Offspring counts are distributed as equally as possible across maternal parents
## to reach the target total (project rule: ≈100 individuals per generation).
##
## @param result    Output from optimise_grid_placement()
## @param pop       Population list
## @param target_n  Target offspring total (default 100)
##
## @return data.frame with columns:
##   maternal_ID, paternal_ID, n_offspring

grid_crosses_from_solution <- function(result, pop, target_n = 100) {
  assignment   <- result$solution
  grid         <- result$grid
  spatial_dist <- build_spatial_distance(grid)
  
  indiv_ids <- as.character(pop$indiv_ID)
  n_indiv   <- length(indiv_ids)
  
  ## Map indiv_ID → location_name
  indiv_to_loc <- setNames(names(assignment), as.character(assignment))
  
  ## Equal offspring per maternal parent, distribute remainder to first N parents
  n_each  <- floor(target_n / n_indiv)
  n_extra <- target_n - n_each * n_indiv
  
  crosses <- vector("list", n_indiv)
  
  for (i in seq_len(n_indiv)) {
    mat_id  <- indiv_ids[i]
    mat_loc <- indiv_to_loc[mat_id]
    
    ## Inverse-distance weights to all other locations
    dists   <- spatial_dist[mat_loc, names(assignment)]
    weights <- ifelse(dists > 0, 1 / dists, 0)
    
    ## Zero out self
    self_loc          <- indiv_to_loc[mat_id]
    weights[self_loc] <- 0
    
    if (sum(weights) == 0) weights <- rep(1 / (n_indiv - 1), n_indiv)
    
    ## Sample one paternal partner weighted by proximity
    pat_loc <- sample(names(weights), size = 1L, prob = weights)
    pat_id  <- assignment[pat_loc]
    
    n_off <- n_each + if (i <= n_extra) 1L else 0L
    
    crosses[[i]] <- data.frame(
      maternal_ID = mat_id,
      paternal_ID = as.character(pat_id),
      n_offspring = n_off,
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, crosses)
}

## ── Apply Mendelian inheritance to produce offspring genotypes ─────────────────
## Given a crossing plan (from grid_crosses_from_solution or scenario scripts),
## generates offspring genotypes using Mendelian inheritance:
##   For each locus, the offspring receives one allele from each parent,
##   sampled from their dosage (0 = AA, 1 = Aa, 2 = aa):
##     dosage → alleles drawn as Binomial(1, dosage/2) from each parent
##
## @param crossing_plan  data.frame with maternal_ID, paternal_ID, n_offspring
## @param pop            Population list supplying parental genotypes ($genotypes)
## @param generation     Character label for the offspring generation (e.g. "F1")
## @param next_id_start  Integer: starting value for new individual IDs
##
## @return A new population list in the standard project format

generate_offspring <- function(crossing_plan, pop, generation = "F1",
                               next_id_start = 1L, scenario_id = NULL) {
  
  ## Handle output directory creation for the specific scenario
  out_dir_path <- "data"
  if (!is.null(scenario_id)) {
    out_dir_path <- paste0("data/scenario_", scenario_id)
    if (!dir.exists(out_dir_path)) {
      dir.create(out_dir_path, recursive = TRUE)
    }
  }
  
  ## Create a tracking log for the crosses
  cross_log_file <- file.path(out_dir_path, paste0("cross_log_", generation, ".csv"))
  
  ## We must handle appending if the file exists (to support tracking across scenarios)
  if (file.exists(cross_log_file)) {
    write.table(crossing_plan, cross_log_file, row.names = FALSE, col.names = FALSE, sep=",", append=TRUE)
  } else {
    write.csv(crossing_plan, cross_log_file, row.names = FALSE)
  }
  
  cat(sprintf("  Logged %d parental crosses for %s to %s\n", nrow(crossing_plan), generation, cross_log_file))
  
  ## Index parental genotypes by ID
  parent_geno <- setNames(pop$genotypes, as.character(pop$indiv_ID))
  n_loci      <- length(pop$genotypes[[1]])
  
  offspring_list <- vector("list", sum(crossing_plan$n_offspring))
  maternal_ids   <- integer(0)
  paternal_ids   <- integer(0)
  count          <- 0L
  
  # Tracking individual IDs along with parental GEBVs
  offspring_id_record <- data.frame(
    Offspring_ID = integer(),
    Maternal_ID = character(),
    Maternal_GEBV = numeric(),
    Paternal_ID = character(),
    Paternal_GEBV = numeric(),
    Generation = character(),
    stringsAsFactors = FALSE
  )
  
  ## Lookup parental GEBVs
  parent_gebv <- setNames(pop$GEBV, as.character(pop$indiv_ID))
  
  
  for (r in seq_len(nrow(crossing_plan))) {
    mat_id  <- as.character(crossing_plan$maternal_ID[r])
    pat_id  <- as.character(crossing_plan$paternal_ID[r])
    n_off   <- crossing_plan$n_offspring[r]
    
    mat_geno <- parent_geno[[mat_id]]
    pat_geno <- parent_geno[[pat_id]]
    
    for (k in seq_len(n_off)) {
      count <- count + 1L
      current_new_id <- next_id_start + count - 1L
      
      ## Mendelian sampling at each locus:
      ## parent dosage d → prob of passing ALT allele = d/2
      mat_allele <- rbinom(n_loci, size = 1L, prob = mat_geno / 2)
      pat_allele <- rbinom(n_loci, size = 1L, prob = pat_geno / 2)
      off_dosage <- mat_allele + pat_allele     # 0, 1, or 2
      
      offspring_list[[count]] <- off_dosage
      maternal_ids            <- c(maternal_ids, as.integer(mat_id))
      paternal_ids            <- c(paternal_ids, as.integer(pat_id))
      
      # Log individual seedling tracking
      offspring_id_record <- rbind(offspring_id_record, data.frame(
        Offspring_ID = current_new_id,
        Maternal_ID = mat_id,
        Maternal_GEBV = parent_gebv[[mat_id]],
        Paternal_ID = pat_id,
        Paternal_GEBV = parent_gebv[[pat_id]],
        Generation = generation,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  n_offspring  <- count
  new_ids      <- seq(next_id_start, length.out = n_offspring)
  
  ## Write out the individual seedling tracking log
  seedling_log_file <- file.path(out_dir_path, paste0("seedling_parent_log_", generation, ".csv"))
  if (file.exists(seedling_log_file)) {
    write.table(offspring_id_record, seedling_log_file, row.names = FALSE, col.names = FALSE, sep=",", append=TRUE)
  } else {
    write.csv(offspring_id_record, seedling_log_file, row.names = FALSE)
  }
  cat(sprintf("  Logged %d individual seedlings and their parents to %s\n", nrow(offspring_id_record), seedling_log_file))
  
  
  list(
    indiv_ID    = new_ids,
    genotypes   = offspring_list,
    GEBV        = rep(NA_real_, n_offspring),   # filled after GAPIT run
    allele_freqs = pop$allele_freqs,            # inherited from original pop
    maternal_ID = maternal_ids,
    paternal_ID = paternal_ids,
    generation  = generation,
    error       = rep(0, n_offspring)
  )
}

## ═══════════════════════════════════════════════════════════════════════════════
## VISUALISATION
## ═══════════════════════════════════════════════════════════════════════════════

## ── Print grid layout ─────────────────────────────────────────────────────────
## Prints the grid with individual IDs at their assigned positions.
## If pop is supplied, also shows GEBV in parentheses.
##
## @param result  Output from optimise_grid_placement()
## @param pop     Optional population list (adds GEBV annotation)

print_grid_layout <- function(result, pop = NULL) {
  grid       <- result$grid
  assignment <- result$solution
  x_max      <- max(grid$x_pos)
  y_max      <- max(grid$y_pos)
  
  ## Build annotation strings
  cell_label <- function(id) {
    if (!is.null(pop)) {
      gebv_val <- pop$GEBV[match(id, as.character(pop$indiv_ID))]
      sprintf("%s(%.3f)", id, gebv_val)
    } else {
      sprintf("%4s", id)
    }
  }
  
  ## Fill a character matrix
  col_width <- if (!is.null(pop)) 12L else 6L
  mat       <- matrix(strrep(" ", col_width), nrow = y_max, ncol = x_max)
  
  for (i in seq_len(nrow(grid))) {
    loc <- grid$location_name[i]
    id  <- assignment[loc]
    x   <- grid$x_pos[i]
    y   <- grid$y_pos[i]
    mat[y, x] <- formatC(cell_label(id), width = col_width, flag = "-")
  }
  
  cat(sprintf("Grid layout [%d × %d] — objective = %.6f\n",
              x_max, y_max, result$final_score))
  cat(strrep("-", col_width * x_max + x_max - 1L), "\n")
  for (y in seq_len(y_max)) {
    cat(paste(mat[y, ], collapse = "|"), "\n")
  }
  cat(strrep("-", col_width * x_max + x_max - 1L), "\n")
  invisible(mat)
}

## ── Plot SA convergence chain ─────────────────────────────────────────────────
## Produces a line plot of score_current vs step number.
## Requires no packages beyond base graphics.
##
## @param result  Output from optimise_grid_placement()

plot_grid_chain <- function(result) {
  chain <- result$chain
  plot(
    chain$step, chain$score_current,
    type  = "l", col  = "steelblue", lwd = 1.5,
    xlab  = "SA step",
    ylab  = "Objective score (lower = better)",
    main  = sprintf(
      "Grid SA convergence (gebv_weight = %.2f, final = %.6f)",
      result$gebv_weight, result$final_score
    )
  )
  lines(chain$step, chain$score_proposed, col = "grey70", lwd = 0.8)
  legend("topright",
         legend = c("Accepted", "Proposed"),
         col    = c("steelblue", "grey70"),
         lwd    = c(1.5, 0.8), bty = "n")
}
