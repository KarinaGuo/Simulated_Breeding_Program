## Scenario 5
## Top 1/3 maternal x Random paternal from Top 1/3 F0 (F1) / all F1 (F2)

source("code/scenario_utils.R")

target_generations <- 5 

stats_log <- list()
stats_log <- append_stats(stats_log, calc_summary_stats(original_init_population, original_init_population))

current_pop <- original_init_population
id_tracker <- length(original_init_population$indiv_ID) + 1

for (gen_num in 1:target_generations) {
  gen_label <- paste0("F", gen_num)
  cat(sprintf("Running generation %s...\n", gen_label))
  
  is_fixed_gen <- (gen_num == 1)
  top_third_pop <- get_top_fraction(current_pop, 1/3, force_fixed_maternal = is_fixed_gen)
  
  if (gen_num == 1) {
    paternal_pop <- top_third_pop
  } else {
    paternal_pop <- current_pop
  }
  
  if (gen_num < target_generations) {
    plan <- random_crossing_plan(top_third_pop, top_third_pop, target_n = 100)
  } else {
    plan <- random_crossing_plan(top_third_pop, paternal_pop, target_n = 100)
  }
  
  if (is_fixed_gen && !is.null(fixed_maternal_pop)) {
    # Combine maternal selection (which now includes the fixed IDs) and the paternal candidates
    var_mat_pop <- top_third_pop
    
    combined_genotypes <- c(var_mat_pop$genotypes, current_pop$genotypes)
    combined_ids <- c(var_mat_pop$indiv_ID, current_pop$indiv_ID)
    combined_gebvs <- c(var_mat_pop$GEBV, current_pop$GEBV)
    
    # Remove duplicates if any
    unique_idx <- !duplicated(combined_ids)
    
    combined_pop_full <- list(
      indiv_ID = combined_ids[unique_idx],
      genotypes = combined_genotypes[unique_idx],
      GEBV = combined_gebvs[unique_idx],
      allele_freqs = original_init_population$allele_freqs
    )
    next_pop <- generate_offspring(plan, combined_pop_full, generation = gen_label, next_id_start = id_tracker, scenario_id = 5)
  } else {
    next_pop <- generate_offspring(plan, current_pop, generation = gen_label, next_id_start = id_tracker, scenario_id = 5)
  }  
  next_pop$GEBV <- Phenotype_from_genotype_GAPIT(next_pop$genotypes, SNPs_tested)
  next_pop$GEBV <- pmax(pmin(next_pop$GEBV, 1), 0.001)
  id_tracker <- id_tracker + length(next_pop$indiv_ID)
  
  stats_log <- append_stats(stats_log, calc_summary_stats(next_pop, original_init_population))
  current_pop <- next_pop
}

stats_df <- stats_to_df(stats_log)
write.csv(stats_df, "data/scenario_5_stats.csv", row.names = FALSE)
saveRDS(current_pop, "data/scenario_5_final_pop.rds")
print(stats_df)