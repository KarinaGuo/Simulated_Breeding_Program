

reconcile_sample_nondominated_archive <- function(summary_1, summary_2, weights, archive) {

   new_archive <- archive

   alen <- length(archive)

   new_dominated   <- FALSE
   i_old_dominated <- c()

   for (i in 1:alen) {

      new_v1_larger <- FALSE; new_v2_larger <- FALSE
      if (summary_1 > archive[[ i ]]$value_1) {new_v1_larger <- TRUE}
      if (summary_2 > archive[[ i ]]$value_2) {new_v2_larger <- TRUE}

      if ( new_v1_larger & new_v2_larger ) {  i_old_dominated <- c(i_old_dominated, i) }

      if ( !new_v1_larger & !new_v2_larger ) {  new_dominated <- TRUE }

   }

   if (!new_dominated) {

      if ( !is.null(i_old_dominated) ) {
         new_archive <- new_archive[ -i_old_dominated ]
      }

      lna <- length(new_archive)
      il  <- lna + 1
      new_archive[[ il ]] <- list(value_1=summary_1, value_2=summary_2, weights=weights)

   }

   return(new_archive)
}


