## pool_results.R
## Pools all scenario results, generates a combined summary table,
## and creates a scatter plot of Mean GEBV and Allelic Capture across generations.


library(tidyverse)

# 1. Load and pool all scenario statistics
all_stats <- list()
for (i in 1:7) {
  file_path <- sprintf("data/scenario_%d_stats.csv", i)
  if (file.exists(file_path)) {
    df <- read.csv(file_path)
    df$Scenario <- paste("Scenario", i)
    all_stats[[i]] <- df
  } else {
    warning(sprintf("File %s not found. Skipping...", file_path))
  }
}

combined_stats <- bind_rows(all_stats)

# Extract numeric generation number for better continuous plotting (F0 = 0, F1 = 1, etc.)
combined_stats$gen_num <- as.numeric(gsub("F", "", combined_stats$generation))

# 2. Generate and save the summary statistic table
cat("=== Combined Summary Statistics Table ===\n")
print(combined_stats)
write.csv(combined_stats, "output/combined_scenario_stats.csv", row.names = FALSE)
cat("\nSaved combined table to output/combined_scenario_stats.csv\n")

# 3. Generate scatter plot with generation on the x-axis, and GEBV / Allelic capture on the y-axis
plot_data <- combined_stats %>%
  pivot_longer(
    cols = c(mean_GEBV, allelic_capture),
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = recode(Metric, 
                    mean_GEBV = "Mean GEBV", 
                    allelic_capture = "Allelic Capture")
  )

# Create the plot
p <- ggplot(plot_data, aes(x = gen_num, y = Value, color = Scenario, group = Scenario)) +
  geom_point() +
  geom_line(alpha = 0.7) +
  facet_wrap(~ Metric, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Generation",
    y = "Score",
    color = "Scenario"
  ) +
  scale_x_continuous(breaks = unique(plot_data$gen_num), labels = paste0("F", unique(plot_data$gen_num)))

# Save the plot
ggsave("output/scenario_comparison_plot.png", plot = p, units = 'px', width = 2000, height = 1500, limitsize = FALSE)
cat("Saved scatter plot to output/scenario_comparison_plot.png\n")
