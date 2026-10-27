## genotype_phenotype_v1.R
##
## Converts a list of individual dosage genotypes to GEBV scores via GAPIT gBLUP.
##
## Requires in the calling environment:
##   hapmap_meta         - first 11 cols of the training hapmap data.frame
##   training_hapmap     - full training-individual hapmap (rs. + S_* cols)
##   ref_alleles         - character vector of REF alleles, one per SNP (ordered)
##   alt_alleles         - character vector of ALT alleles, one per SNP (ordered)
##   het_codes           - named character vector of IUPAC heterozygote codes
##   input_phenodatafile - path to phenotype CSV
##   out_dir             - output directory for temp hapmap and GAPIT predictions
##
## @param individuals_GT  List of N integer vectors (dosage 0/1/2 per SNP).
##                        This is `population$genotypes`.
## @param SNPs_tested     Character vector of SNP names (must match hapmap_meta$rs.)
##
## @return Numeric vector of length N, GEBV rescaled to [0.001, 1].
##         Lower values = more resistant (better).

Phenotype_from_genotype_GAPIT <- function(individuals_GT, SNPs_tested) {

  ## ── Build dosage matrix ────────────────────────────────────────────────────
  gt_matrix  <- do.call(cbind, individuals_GT)      # n_SNPs x N
  sim_names  <- paste0("Sim_", seq_len(ncol(gt_matrix)))
  colnames(gt_matrix) <- sim_names

  ## ── Encode dosage → IUPAC ─────────────────────────────────────────────────
  iupac_matrix <- matrix("N", nrow = nrow(gt_matrix), ncol = ncol(gt_matrix))
  colnames(iupac_matrix) <- sim_names

  iupac_matrix[gt_matrix == 0] <- ref_alleles[row(gt_matrix)[gt_matrix == 0]]
  iupac_matrix[gt_matrix == 1] <- het_codes[row(gt_matrix)[gt_matrix == 1]]
  iupac_matrix[gt_matrix == 2] <- alt_alleles[row(gt_matrix)[gt_matrix == 2]]
  iupac_matrix[is.na(iupac_matrix)] <- "N"

  iupac_df <- data.frame(
    rs. = SNPs_tested,
    iupac_matrix,
    stringsAsFactors = FALSE
  )

  ## ── Join with hapmap metadata ──────────────────────────────────────────────
  iupac_df$rs.      <- as.character(iupac_df$rs.)
  hapmap_meta$rs.   <- as.character(hapmap_meta$rs.)
  simulated_hapmap  <- left_join(iupac_df, hapmap_meta, by = "rs.")

  ## ── Concatenate with training individuals ──────────────────────────────────
  simulated_hapmap_trainingconcat <- left_join(training_hapmap, simulated_hapmap, by = "rs.")

  ## ── Write temp hapmap ─────────────────────────────────────────────────────
  tmp_hapmap_path <- file.path(out_dir, "tmp_simulated.hapmap.hmp.txt")
  write.table(simulated_hapmap_trainingconcat, file = tmp_hapmap_path,
              row.names = FALSE, sep = "\t")
  gt_datafile_tmp <- normalizePath(tmp_hapmap_path, mustWork = FALSE)

  ## ── Run GAPIT ─────────────────────────────────────────────────────────────
  filename <- "/prediction_gt_simulated_run.csv"

  run_genomic_prediction_simulation(
    gt_datafile         = gt_datafile_tmp,
    input_phenodatafile = input_phenodatafile,
    trans_COI           = TRUE,
    rm_clim_snp         = TRUE,
    homo_site_filt      = 1,
    miss_site_filt      = 1,
    homo_samp_filt      = 1,
    miss_samp_filt      = 1,
    out_dir             = out_dir,
    filename            = filename
  )

  ## ── Extract and rescale predictions ───────────────────────────────────────
  pred     <- read.csv(file.path(out_dir, "prediction_gt_simulated_run.csv"))
  pred_sim <- pred[pred$Taxa %in% sim_names, ]
  pred_sim$Taxa <- factor(pred_sim$Taxa, levels = sim_names)
  pred_sim <- pred_sim[order(pred_sim$Taxa), ]

  ## Rescale from raw GAPIT prediction range [0,3] to trait range [0.001, 1]
  ## (lower = more resistant = lower COI)
  MR <- scales::rescale(pred_sim$Prediction, from = c(0, 3), to = c(0.001, 1))

  return(MR)
}
