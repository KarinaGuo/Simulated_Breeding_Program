calc_QC_gt <- function(gt_data=myG, homo_site_filt, miss_site_filt, miss_samp_filt, homo_samp_filt, rm_clim_snp, log_file){
  log_appender(appender_file(file = log_file))
  myG <- gt_data
  
  ## Generating site and indiv filtering data wrangle
  genetics_data_filtering <- myG
  colnames(genetics_data_filtering) <- genetics_data_filtering[1,]; genetics_data_filtering <- genetics_data_filtering[c(-1),]; genetics_data_filtering <- genetics_data_filtering[,-c(2:11)]
  genetics_data_t <- t(genetics_data_filtering)
  colnames(genetics_data_t) <- genetics_data_t[1,]; genetics_data_t <- genetics_data_t[c(-1),]
  genetics_data_t <- cbind(sample=row.names(genetics_data_t), genetics_data_t)
  
  genetics_data_t_long <- as.data.frame(genetics_data_t) %>% 
    pivot_longer(
      cols = -sample,   # Keep 'sample' column fixed
      names_to = "SNP", # New column for SNP names
      values_to = "Genotype" # New column for values
    )
  
  genetics_data_t_long_rmWGS <- genetics_data_t_long[!grepl("^S_(\\d+)_.*$", genetics_data_t_long$sample),] # Which individuals are from WGS
  numsamps <- length(unique(genetics_data_t_long_rmWGS$sample))
  numloci <- length(unique(genetics_data_t_long_rmWGS$SNP))
  
  
  unique(genetics_data_t_long_rmWGS$Genotype) # How many unique genotypes, see hapmap encoding above
  
  ### By hapmap breakdown into hetero and homo
  genetics_data_t_long_rmWGS$Zygosity=0 #init
  genetics_data_t_long_rmWGS$Zygosity <- ifelse(genetics_data_t_long_rmWGS$Genotype %in% c("A", "T", "G", "C"), "homozygous", "heterozygous")
  genetics_data_t_long_rmWGS$Zygosity[genetics_data_t_long_rmWGS$Genotype=="N"] <- "missing"
  
  ################
  ### Count incidence
  homo_hetero_samp <- genetics_data_t_long_rmWGS %>%
    group_by(sample) %>% 
    summarise(Total=n(), homo=sum(Zygosity=="homozygous"), hetero=sum(Zygosity=="heterozygous"), miss=sum(Zygosity=="missing")) %>%
    ungroup()
  
  homo_hetero_samp$homo_props <- homo_hetero_samp$homo/(homo_hetero_samp$homo+homo_hetero_samp$hetero)
  homo_hetero_samp$hetero_props <- homo_hetero_samp$hetero/(homo_hetero_samp$homo+homo_hetero_samp$hetero)
  homo_hetero_samp$miss_props <- homo_hetero_samp$miss/homo_hetero_samp$Total
  
  #
  
  homo_hetero_loci <- genetics_data_t_long_rmWGS %>%
    group_by(SNP) %>% 
    summarise(Total=n(), homo=sum(Zygosity=="homozygous"), hetero=sum(Zygosity=="heterozygous"), miss=sum(Zygosity=="missing")) %>%
    ungroup()
  
  homo_hetero_loci$homo_props <- homo_hetero_loci$homo/(homo_hetero_loci$homo + homo_hetero_loci$hetero)
  homo_hetero_loci$hetero_props <- homo_hetero_loci$hetero/(homo_hetero_loci$homo + homo_hetero_loci$hetero)
  homo_hetero_loci$miss_props <- homo_hetero_loci$miss/homo_hetero_loci$Total
  
  rm_samps_miss <- (homo_hetero_samp %>% filter (miss_props > miss_samp_filt))
  log_info("Removing high miss samps ", nrow(rm_samps_miss), " out of ", numsamps, "\n")
  
  rm_loci_miss <- (homo_hetero_loci %>% filter (miss_props >  miss_site_filt))
  log_info("Removing high miss loci ", nrow(rm_loci_miss), " out of ", numloci, "\n")
  
  rm_samps_homo <- (homo_hetero_samp %>% filter (homo_props > homo_samp_filt))
  log_info("Removing high homo samps ", nrow(rm_samps_homo), " out of ", numsamps, "\n")
  
  rm_loci_homo <- (homo_hetero_loci %>% filter (homo_props > homo_site_filt))
  log_info("Removing high homo loci ", nrow(rm_loci_homo), " out of ", numloci, "\n")
  
  if(rm_clim_snp == TRUE){rm_loci_MR <- (homo_hetero_loci[grepl("_climate", homo_hetero_loci$SNP),])}
  
  rm_samps <- unique(c(rm_samps_homo$sample, rm_samps_miss$sample))
  rm_loci <- unique(append(rm_loci_homo$SNP, rm_loci_miss$SNP)); 
  
  if(rm_clim_snp == TRUE){rm_loci <- append(rm_loci,rm_loci_MR$SNP)}
  
  log_info("Total samps removed: ", length(rm_samps), "\n")
  log_info("Total loci removed: ", length(rm_loci), "\n")
  
  ### Removing samps and loci
  myG_rmloci <- myG %>% dplyr::filter(!(V1 %in% rm_loci))
  myG_rmloci_samps <- myG_rmloci[, !(myG_rmloci[1,] %in% rm_samps)]
  
  return(myG_rmloci_samps)
}

run_genomic_prediction <- function (gt_datafile, input_phenodatafile, remove_SNPs=FALSE, SNPS_TBR=NULL, trans_COI=TRUE, rm_clim_snp=FALSE, out_dir=getwd(), homo_site_filt=0.8, miss_site_filt=0.9, homo_samp_filt=0.9, miss_samp_filt=0.9, log_file, filename ="prediction_gt.csv"){
  require(logger)
  log_appender(appender_file(file = log_file))
  
  log_info("----------------------- New run ----------------------------- \n")
  log_info("Output directory set to: ",out_dir,"\n")
  log_info("Prediction ouput file set to: ",paste0(out_dir,"/prediction_gt.csv \n"))
  
  log_info("Running with filters of: ", 
           " \n Individual missingness: ", miss_samp_filt, 
           " \n Individual homozygosity: ", homo_samp_filt, 
           " \n Site missingness: ", miss_site_filt, 
           " \n Site homozygosity: ", homo_site_filt, "\n")
  
  require(tidyverse)
  require(patchwork)
  
  # tryCatch({
  #  source("http://zzlab.net/GAPIT/gapit_functions.txt")
  #}, error = function(e) {
  #  message("Original source failed, using local backup...")
    source("~/Uni/Doctorate/gapit_functions_080425.txt")
  #})
  
  #Import files
  seedling_mquin_rust <- read.csv(input_phenodatafile, sep = ",") # Rust data of seedling individuals
  colnames(seedling_mquin_rust) <- c('LIBRARY', 'COI')
  
  Y <- seedling_mquin_rust
  if(trans_COI==TRUE){Y$COI <- (Y$COI)^0.25}
  # Replacing Chrom names
  chromosome_mapping <- c("MQA_CHR01" = 1, "MQA_CHR02" = 2,"MQA_CHR03" = 3, "MQA_CHR04" = 4, "MQA_CHR05" = 5, "MQA_CHR06" = 6, "MQA_CHR07" = 7, "MQA_CHR08" = 8, "MQA_CHR09" = 9, "MQA_CHR10" = 10, "MQA_CHR11" = 11)
  
  file_dir = out_dir
  if(!file.exists(file_dir)){
    print(paste("Directory", file_dir, "does not exist and is being created :)")) 
    dir.create(file_dir)
  }
  wd = file_dir
  
  # Running gBLUP for GP
  setwd(wd)
  myG <- read.csv(file = gt_datafile, sep = "\t", header = FALSE)
  myG$V3 <- chromosome_mapping[myG$V3] 
  myG$V3[1] <- "chrom"
  
  myG <- calc_QC_gt(gt_data=myG, homo_site_filt=homo_site_filt, miss_site_filt=miss_site_filt, miss_samp_filt=miss_samp_filt, homo_samp_filt=homo_samp_filt, log_file=log_file, rm_clim_snp=rm_clim_snp)
  
  if (remove_SNPs) {
    myG <- myG %>% dplyr::filter(!(V1 %in% SNPS_TBR))
  }
  
  proceed <- "Y"
  
  if (tolower(trimws(proceed)) == "y") {
    myGAPIT <- tryCatch(expr = {(GAPIT(Y=Y, G=myG, PCA.total=3, model=c("gBLUP"), file.output=FALSE))}, error = function(e){print("Running Geno.View.output = FALSE"); (GAPIT(Y=Y, G=myG, PCA.total=3, model=c("gBLUP"), Geno.View.output = FALSE))})
    
    # Extracting predictions + plotting
    prediction=myGAPIT$Pred
    unique(prediction$RefInf)
    
    colnames(Y)[1] <- "Taxa"
    prediction_gt <- left_join(prediction, Y)
    prediction_gt$Prediction[prediction_gt$Prediction < 0] <- 0
    
    log_info("GAPIT h2 val=", myGAPIT$h2)
    
    cat("Saving file to:", paste0(wd,filename))
    write.csv(prediction_gt, file = paste0(wd,filename))
  } else {
    cat("Run cancelled.\n"); log_info("Run cancelled.\n")
    return(NULL)
  }
}

calc_QC_gt <- function(gt_data=myG, homo_site_filt, miss_site_filt, miss_samp_filt, homo_samp_filt, rm_clim_snp){
  log_appender(appender_file(file = log_file))
  myG <- gt_data
  
  ## Generating site and indiv filtering data wrangle
  genetics_data_filtering <- myG
  colnames(genetics_data_filtering) <- genetics_data_filtering[1,]; genetics_data_filtering <- genetics_data_filtering[c(-1),]; genetics_data_filtering <- genetics_data_filtering[,-c(2:11)]
  genetics_data_t <- t(genetics_data_filtering)
  colnames(genetics_data_t) <- genetics_data_t[1,]; genetics_data_t <- genetics_data_t[c(-1),]
  genetics_data_t <- cbind(sample=row.names(genetics_data_t), genetics_data_t)
  
  genetics_data_t_long <- as.data.frame(genetics_data_t) %>% 
    pivot_longer(
      cols = -sample,   # Keep 'sample' column fixed
      names_to = "SNP", # New column for SNP names
      values_to = "Genotype" # New column for values
    )
  
  genetics_data_t_long_rmWGS <- genetics_data_t_long[!grepl("^S_(\\d+)_.*$", genetics_data_t_long$sample),] # Which individuals are from WGS
  numsamps <- length(unique(genetics_data_t_long_rmWGS$sample))
  numloci <- length(unique(genetics_data_t_long_rmWGS$SNP))
  
  
  unique(genetics_data_t_long_rmWGS$Genotype) # How many unique genotypes, see hapmap encoding above
  
  ### By hapmap breakdown into hetero and homo
  genetics_data_t_long_rmWGS$Zygosity=0 #init
  genetics_data_t_long_rmWGS$Zygosity <- ifelse(genetics_data_t_long_rmWGS$Genotype %in% c("A", "T", "G", "C"), "homozygous", "heterozygous")
  genetics_data_t_long_rmWGS$Zygosity[genetics_data_t_long_rmWGS$Genotype=="N"] <- "missing"
  
  ################
  ### Count incidence
  homo_hetero_samp <- genetics_data_t_long_rmWGS %>%
    group_by(sample) %>% 
    summarise(Total=n(), homo=sum(Zygosity=="homozygous"), hetero=sum(Zygosity=="heterozygous"), miss=sum(Zygosity=="missing")) %>%
    ungroup()
  
  homo_hetero_samp$homo_props <- homo_hetero_samp$homo/(homo_hetero_samp$homo+homo_hetero_samp$hetero)
  homo_hetero_samp$hetero_props <- homo_hetero_samp$hetero/(homo_hetero_samp$homo+homo_hetero_samp$hetero)
  homo_hetero_samp$miss_props <- homo_hetero_samp$miss/homo_hetero_samp$Total
  
  #
  
  homo_hetero_loci <- genetics_data_t_long_rmWGS %>%
    group_by(SNP) %>% 
    summarise(Total=n(), homo=sum(Zygosity=="homozygous"), hetero=sum(Zygosity=="heterozygous"), miss=sum(Zygosity=="missing")) %>%
    ungroup()
  
  homo_hetero_loci$homo_props <- homo_hetero_loci$homo/(homo_hetero_loci$homo + homo_hetero_loci$hetero)
  homo_hetero_loci$hetero_props <- homo_hetero_loci$hetero/(homo_hetero_loci$homo + homo_hetero_loci$hetero)
  homo_hetero_loci$miss_props <- homo_hetero_loci$miss/homo_hetero_loci$Total
  
  rm_samps_miss <- (homo_hetero_samp %>% filter (miss_props > miss_samp_filt))
  log_info("Removing high miss samps ", nrow(rm_samps_miss), " out of ", numsamps, "\n")
  
  rm_loci_miss <- (homo_hetero_loci %>% filter (miss_props >  miss_site_filt))
  log_info("Removing high miss loci ", nrow(rm_loci_miss), " out of ", numloci, "\n")
  
  rm_samps_homo <- (homo_hetero_samp %>% filter (homo_props > homo_samp_filt))
  log_info("Removing high homo samps ", nrow(rm_samps_homo), " out of ", numsamps, "\n")
  
  rm_loci_homo <- (homo_hetero_loci %>% filter (homo_props > homo_site_filt))
  log_info("Removing high homo loci ", nrow(rm_loci_homo), " out of ", numloci, "\n")
  
  if(rm_clim_snp == TRUE){rm_loci_MR <- (homo_hetero_loci[grepl("_climate", homo_hetero_loci$SNP),])}
  
  rm_samps <- unique(c(rm_samps_homo$sample, rm_samps_miss$sample))
  rm_loci <- unique(append(rm_loci_homo$SNP, rm_loci_miss$SNP)); 
  
  if(rm_clim_snp == TRUE){rm_loci <- append(rm_loci,rm_loci_MR$SNP)}
  
  log_info("Total samps removed: ", length(rm_samps), "\n")
  log_info("Total loci removed: ", length(rm_loci), "\n")
  
  ### Removing samps and loci
  myG_rmloci <- myG %>% dplyr::filter(!(V1 %in% rm_loci))
  myG_rmloci_samps <- myG_rmloci[, !(myG_rmloci[1,] %in% rm_samps)]
  
  return(myG_rmloci_samps)
}

run_genomic_prediction_simulation <- function (gt_datafile, input_phenodatafile, remove_SNPs=FALSE, SNPS_TBR=NULL, trans_COI=TRUE, rm_clim_snp=FALSE, out_dir=getwd(), homo_site_filt=0.8, miss_site_filt=0.9, homo_samp_filt=0.9, miss_samp_filt=0.9, log_file, filename ="prediction_gt.csv"){
  #Import files
  seedling_mquin_rust <- read.csv(input_phenodatafile, sep = ",") # Rust data of seedling individuals
  colnames(seedling_mquin_rust) <- c('LIBRARY', 'COI')
  
  Y <- seedling_mquin_rust
  if(trans_COI==TRUE){Y$COI <- (Y$COI)^0.25}
  # Replacing Chrom names
  chromosome_mapping <- c("MQA_CHR01" = 1, "MQA_CHR02" = 2,"MQA_CHR03" = 3, "MQA_CHR04" = 4, "MQA_CHR05" = 5, "MQA_CHR06" = 6, "MQA_CHR07" = 7, "MQA_CHR08" = 8, "MQA_CHR09" = 9, "MQA_CHR10" = 10, "MQA_CHR11" = 11)
  
  file_dir = out_dir
  if(!file.exists(file_dir)){
    print(paste("Directory", file_dir, "does not exist and is being created :)")) 
    dir.create(file_dir)
  }
  og_wd <- getwd()
  wd = file_dir
  
  # Running gBLUP for GP
  setwd(wd)
  myG <- read.csv(file = gt_datafile, sep = "\t", header = FALSE)
  chrom_idx = match("chrom", myG[1,])
  myG[,chrom_idx] <- chromosome_mapping[myG[,chrom_idx]] 
  myG[1, chrom_idx] <- "chrom"
  
  myG <- myG
  
  if (remove_SNPs) {
    myG <- myG %>% dplyr::filter(!(V1 %in% SNPS_TBR))
  }
  
  geno_calls <- myG[2:nrow(myG), 12:ncol(myG)]
  is_clone <- duplicated(t(geno_calls))
  clone_cols <- which(is_clone) + 11
  for (col in clone_cols) {
    current_call <- myG[2, col]
    myG[2, col] <- ifelse(current_call == "A", "T", "A")
  }
  
  
  proceed <- "Y"
  
  if (tolower(trimws(proceed)) == "y") {
    myGAPIT <- tryCatch(expr = {(GAPIT(Y=Y, G=myG, model=c("gBLUP"),file.output=FALSE))}, error = function(e){print("Running Geno.View.output = FALSE"); (GAPIT(Y=Y, G=myG, model=c("gBLUP"), Geno.View.output = FALSE))})
    
    # Extracting predictions + plotting
    prediction=myGAPIT$Pred
    
    colnames(Y)[1] <- "Taxa"
    prediction_gt <- left_join(prediction, Y)
    #prediction_gt$Prediction[prediction_gt$Prediction < 0] <- 0
    
    write.csv(prediction_gt, file = paste0(getwd(),filename))
  } else {
    cat("Run cancelled.\n")
    return(NULL)
  }
  
  setwd(og_wd)
}

  