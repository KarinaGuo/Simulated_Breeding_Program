# Project Overview

## Inputs

* 30 individuals with random genotypes for the snps in SNP_list.csv are generated 

## Rules
* If both a maternal and paternal individual is identified, individuals inherit the genotype using Mendelian rules
* If only a maternal individual is identified, the paternal individual is identified using the rules in the scenario's outline. The individual then inherits the genotype using Mendelian rules
* The trait is limited to 0.001 to 1 
* Each individual has an equal number of offsprings, with the final number of offsprings to be as close as possible to a total of 100 individuals
* At the end of each scenario, track the number of generations and the summary statistic scores.
* Create a new script for each scenario
* The lower the GEBV score is optimal

## Assumptions
* Mendelian inheritance
* No drift/mutation
* No mortality or immigration of new individuals

## Initial things to do
* Outline input and output data structure
* Ensure code\initialise_population.R is operating correctly
* Create scripts to generate summary statistic scores. Scores include, allelic diversity within population, allelic capture from the original_init_population's genotype, mean GEBV, GEBV standard deviation.
* Create a script to calculate kinship from the genotypes of each individual
* Read and refactor grid optimisation script to accommodate for input and outputs of this project (optimisation_scripts.R and optimisation_workflow.R)
* Refactor kinship + GEBV minimiser optimiser to accommodate for input and outputs of this project (OptGenMix-master)


## Implementation plan
* 1. Generate the initial population genotypes (original_init_population) using code\initialise_population.R
* 2. Run Phenotype_from_genotype_GAPIT on init_AF to generate GEBV scores for all individuals
* 3. Select individuals for crossing as outlined in the scenario.
* 4. Seedlings are generated with the genotype from their selected parental lines.
* 5. Run Phenotype_from_genotype_GAPIT to generate GEBV scores for all new individuals
* 6. Generate summary statistic scores
* 7. Cross individuals as outlined in the scenario
* 8. Rerun 4-6
* 9. If the scenario dictates it, repeat for the number of generations indicated

## Scenarios
* Capacity to upload F0 
* Capacity to fix maternal lines to be used in the final generation - fixing is not integrated into Scenarios 7-8
* Backcrossing occurs with the whole population of the set generation

### Scenario 1
* Subset the top 1/3 of the GEBV from the original_init_population for the initial cross population, this is henceforth called the top_third_population. Select each individual in the top_third_population as the maternal line. For each maternal line, select a paternal line randomly from the top_third_population, with repeats allowed. This produces the F1 population.
* Generate the summary statistics of the F1 population
* Repeat the above until F2

### Scenario 2
* Subset the top 1/2 of the GEBV from the original_init_population for the initial cross population. Repeat Scenario 1.

### Scenario 3
* Subset the top 1/3 of the GEBV from the original_init_population for the initial cross population, this is henceforth called the top_third_population. Select each individual in the top_third_population as the maternal line.
* Select the paternal line from the top_third_population by minimising kinship and GEBV against the maternal line

### Scenario 4
* Subset the top 1/2 of the GEBV from the original_init_population for the initial cross population. Repeat Scenario 3.

### Scenario 5-6
* Repeat Scenario 1-2 but select the paternal line from the F1 for the last generation

### Scenario 7
* Cross all individuals minimising only relatedness

### Scenario 8
* Generate a grid that fits all the individuals with the shortest total edge length possible
* Place all individuals within the grid to optimise for the minimal kinship and minimal GEBV across all individuals. 
* Cross individuals weighted by grid distance.