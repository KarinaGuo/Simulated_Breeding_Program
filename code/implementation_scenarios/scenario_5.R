## Scenario 5
## Top 1/3 maternal x Random paternal from Top 1/3 F0 (F1) / all F1 (F2)

source("code/scenario_utils.R")

target_generations <- 5 

stats_log <- list()
stats_log <- append_stats(stats_log, calc_summary_stats(original_init_population, original_init_population))

current_pop <- original_init_population

for (gen_num in 1:target_generations) {
  gen_label <- paste0("F", gen_num)
  cat(sprintf("Running generation %s...\n", gen_label))
  
  top_third_pop <- get_top_fraction(current_pop, 1/3)
  
  if (gen_num == 1) {
    paternal_pop <- top_third_pop
  } else {
    paternal_pop <- current_pop
  }
  
  plan <- random_crossing_plan(top_third_pop, paternal_pop, target_n = 100)
  
  next_id <- max(current_pop$indiv_ID) + 1
  next_pop <- generate_offspring(plan, current_pop, generation = gen_label, next_id_start = next_id)
  
  next_pop$GEBV <- Phenotype_from_genotype_GAPIT(next_pop$genotypes, SNPs_tested)
  next_pop$GEBV <- pmax(pmin(next_pop$GEBV, 1), 0.001)
  
  stats_log <- append_stats(stats_log, calc_summary_stats(next_pop, original_init_population))
  current_pop <- next_pop
}

stats_df <- stats_to_df(stats_log)
write.csv(stats_df, "data/scenario_5_stats.csv", row.names = FALSE)
saveRDS(current_pop, "data/scenario_5_final_pop.rds")
print(stats_df)