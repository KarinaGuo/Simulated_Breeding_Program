library(shiny)
library(tidyverse)
library(ggplot2)

# Load existing backend utilities
setwd("~/Uni/Doctorate/Breeding_project")

source("code/scenario_utils.R")
source("code/summary_statistics.R")
source("code/genotype_phenotype_v1.R")
source("code/GP_functs.R")

# Custom wrapper to handle Fixed Integer selection
get_top_n <- function(pop, n, force_fixed_maternal = FALSE) {
  fraction <- min(1, n / length(pop$indiv_ID))
  get_top_fraction(pop, fraction, force_fixed_maternal)
}

ui <- fluidPage(
  titlePanel("Breeding Project Simulation Interface"),
  
  sidebarLayout(
    sidebarPanel(
      tags$b("1. Genotype Inputs"),
      fileInput("fixed_mat_file", "Upload fixed_maternal_genotypes.csv", accept = ".csv"),
      fileInput("f0_file", "Upload uploaded_F0_genotypes.csv", accept = ".csv"),
      checkboxInput("use_uploaded", "Use Uploaded Files (if available)", value = FALSE),
      
      tags$hr(),
      tags$b("2. Selection Parameters"),
      radioButtons("sel_method", "Selection Method:", 
                   choices = c("Percentage (%)" = "pct", "Fixed Integer (N)" = "int")),
      
      conditionalPanel(
        condition = "input.sel_method == 'pct'",
        numericInput("mat_pct", "Top Maternal Fraction (e.g., 0.33 for 1/3):", value = 0.33, min = 0.01, max = 1, step = 0.01),
        numericInput("pat_pct", "Top Paternal Fraction:", value = 0.33, min = 0.01, max = 1, step = 0.01)
      ),
      
      conditionalPanel(
        condition = "input.sel_method == 'int'",
        numericInput("mat_n", "Top Maternal Integer (N):", value = 10, min = 1, step = 1),
        numericInput("pat_n", "Top Paternal Integer (N):", value = 10, min = 1, step = 1)
      ),
      
      checkboxInput("mat_only", "Only maternal lines selected (Paternal random across whole generation)", value = FALSE),
      
      tags$hr(),
      tags$b("3. Simulation Settings"),
      numericInput("num_seedlings", "Final Number of Seedlings per Generation:", value = 100, min = 1, step = 1),
      numericInput("num_gens", "Generations to Iterate:", value = 5, min = 1, step = 1),
      
      actionButton("run_sim", "Run Simulation", class = "btn-primary", style = "width: 100%; margin-top: 15px;")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Breeding Design Log", verbatimTextOutput("design_log")),
        tabPanel("Simulation Plots", 
                 plotOutput("gebv_plot"),
                 plotOutput("capture_plot")
        ),
        tabPanel("Summary Data", dataTableOutput("stats_table"))
      )
    )
  )
)

server <- function(input, output, session) {
  
  sim_results <- reactiveVal(NULL)
  sim_log <- reactiveVal("Awaiting simulation run...")
  
  observeEvent(input$run_sim, {
    req(input$num_gens, input$num_seedlings)
    
    withProgress(message = 'Simulation Running...', value = 0, {
      
      log_text <- c("Starting simulation run...")
      sim_log(paste(tail(log_text, 5), collapse = "\n"))
      
      # 1. Initialize Baseline Population 
      if (!is.null(input$f0_file) && input$use_uploaded) {
        uploaded_f0_df <<- read.csv(input$f0_file$datapath, check.names = FALSE)
        log_text <- c(log_text, "Using uploaded F0 genotypes.")
      } else {
        uploaded_f0_df <<- NULL
      }
      
      if (!is.null(input$fixed_mat_file) && input$use_uploaded) {
        fixed_maternal_df <<- read.table(input$fixed_mat_file$datapath, sep = ",", header = TRUE, check.names = FALSE)
        log_text <- c(log_text, "Using uploaded fixed maternal genotypes.")
      } else {
        fixed_maternal_df <<- NULL
      }
      sim_log(paste(tail(log_text, 5), collapse = "\n"))
      
      # Source initialization script (assumes it reads uploaded_f0_df and fixed_maternal_df from env)
      source("code/implementation_scenarios/initialise_population.R", local = TRUE)
      
      current_pop <- original_init_population
      id_tracker <- length(original_init_population$indiv_ID) + 1
      stats_log <- list()
      stats_log <- append_stats(stats_log, calc_summary_stats(original_init_population, original_init_population))
      
      # 2. Iterate Generations
      for (gen_num in 1:input$num_gens) {
        gen_label <- paste0("F", gen_num)
        
        is_fixed_gen <- (gen_num == 1)
        
        # Maternal Selection
        if (input$sel_method == "pct") {
          mat_pop <- get_top_fraction(current_pop, input$mat_pct, force_fixed_maternal = is_fixed_gen)
          log_text <- c(log_text, sprintf("Running %s... Selected top %.2f%% for maternal pool.", gen_label, input$mat_pct * 100))
        } else {
          mat_pop <- get_top_n(current_pop, input$mat_n, force_fixed_maternal = is_fixed_gen)
          log_text <- c(log_text, sprintf("Running %s... Selected top %d individuals for maternal pool.", gen_label, input$mat_n))
        }
        
        # Paternal Selection
        if (input$mat_only) {
          pat_pop <- current_pop
          log_text <- c(log_text, "Paternal lines selected randomly across whole generation.")
        } else {
          if (input$sel_method == "pct") {
            pat_pop <- get_top_fraction(current_pop, input$pat_pct, force_fixed_maternal = FALSE)
          } else {
            pat_pop <- get_top_n(current_pop, input$pat_n, force_fixed_maternal = FALSE)
          }
        }
        
        # Update Log (Last 5 lines) and loading widget dynamically
        sim_log(paste(tail(log_text, 5), collapse = "\n"))
        incProgress(1/input$num_gens, detail = paste("Processing", gen_label))
        
        # Crossing and Offspring Generation
        plan <- random_crossing_plan(mat_pop, pat_pop, target_n = input$num_seedlings)
        
        if (is_fixed_gen && !is.null(fixed_maternal_pop)) {
          combined_genotypes <- c(mat_pop$genotypes, current_pop$genotypes)
          combined_ids <- c(mat_pop$indiv_ID, current_pop$indiv_ID)
          combined_gebvs <- c(mat_pop$GEBV, current_pop$GEBV)
          unique_idx <- !duplicated(combined_ids)
          
          combined_pop_full <- list(
            indiv_ID = combined_ids[unique_idx],
            genotypes = combined_genotypes[unique_idx],
            GEBV = combined_gebvs[unique_idx],
            allele_freqs = original_init_population$allele_freqs
          )
          next_pop <- generate_offspring(plan, combined_pop_full, generation = gen_label, next_id_start = id_tracker, scenario_id = 99)
        } else {
          next_pop <- generate_offspring(plan, current_pop, generation = gen_label, next_id_start = id_tracker, scenario_id = 99)
        }  
        
        # Predict GEBV and update trackers
        next_pop$GEBV <- Phenotype_from_genotype_GAPIT(next_pop$genotypes, SNPs_tested)
        next_pop$GEBV <- pmax(pmin(next_pop$GEBV, 1), 0.001)
        id_tracker <- id_tracker + length(next_pop$indiv_ID)
        
        stats_log <- append_stats(stats_log, calc_summary_stats(next_pop, original_init_population))
        current_pop <- next_pop
      }
      
      # 3. Finalize Output Data
      stats_df <- stats_to_df(stats_log)
      stats_df$gen_num <- as.numeric(gsub("F", "", stats_df$generation))
      
      sim_results(stats_df)
      
      # Print full log at completion
      log_text <- c(log_text, "\nSimulation complete.")
      sim_log(paste(log_text, collapse = "\n"))
      
    })
  })
  
  # Outputs
  output$design_log <- renderText({ sim_log() })
  
  output$stats_table <- renderDataTable({ 
    validate(need(sim_results(), "Awaiting simulation run..."))
    sim_results() %>% 
      mutate(across(where(is.numeric), ~ round(., 2)))
  })
  
  output$gebv_plot <- renderPlot({
    validate(need(sim_results(), "Awaiting simulation run..."))
    df <- sim_results()
    ggplot(df, aes(x = gen_num, y = mean_GEBV)) +
      geom_point() +
      geom_line() +
      theme_bw() +
      labs(title = "Mean GEBV Across Generations", x = "Generation", y = "Mean GEBV") +
      scale_x_continuous(breaks = unique(df$gen_num), labels = paste0("F", unique(df$gen_num)))
  })
  
  output$capture_plot <- renderPlot({
    validate(need(sim_results(), "Awaiting simulation run..."))
    df <- sim_results()
    ggplot(df, aes(x = gen_num, y = allelic_capture)) +
      geom_point() +
      geom_line() +
      theme_bw() +
      labs(title = "Allelic Capture Across Generations", x = "Generation", y = "Allelic Capture (Proportion)") +
      scale_x_continuous(breaks = unique(df$gen_num), labels = paste0("F", unique(df$gen_num)))
  })
}

shinyApp(ui, server)