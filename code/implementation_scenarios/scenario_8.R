## Scenario 8
## Grid-based spatial optimisation for crosses

source("code/scenario_utils.R")

target_generations <- 5 

stats_log <- list()
stats_log <- append_stats(stats_log, calc_summary_stats(original_init_population, original_init_population))

current_pop <- original_init_population

for (gen_num in 1:target_generations) {
  gen_label <- paste0("F", gen_num)
  cat(sprintf("Running generation %s...\n", gen_label))
  
  K_mat <- calc_kinship(current_pop, allele_freqs = original_init_population$allele_freqs)
  
  res <- optimise_grid_placement(current_pop, K_mat, gebv_weight = 0.5, max_steps = 50000)
  plan <- grid_crosses_from_solution(res, current_pop, target_n = 100)
  
  next_id <- max(current_pop$indiv_ID) + 1
  next_pop <- generate_offspring(plan, current_pop, generation = gen_label, next_id_start = next_id)
  
  next_pop$GEBV <- Phenotype_from_genotype_GAPIT(next_pop$genotypes, SNPs_tested)
  next_pop$GEBV <- pmax(pmin(next_pop$GEBV, 1), 0.001)
  
  stats_log <- append_stats(stats_log, calc_summary_stats(next_pop, original_init_population))
  current_pop <- next_pop
}

stats_df <- stats_to_df(stats_log)
write.csv(stats_df, "data/scenario_8_stats.csv", row.names = FALSE)
saveRDS(current_pop, "data/scenario_8_final_pop.rds")
print(stats_df)