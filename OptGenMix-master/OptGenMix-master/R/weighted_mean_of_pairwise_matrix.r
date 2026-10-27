#' Calculates the mean value for pairwise matrix, 
#' when supplied with weights that should be 
#' applied to the individuals 
#' 
#' @param sm  - pairwise matrix of values
#' @param w   - vector of weights
#' @return mean similarity 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export

weighted_mean_of_pairwise_matrix <- function(sm, w=NULL) {

   if(!is.matrix(sm)) {
      sm <- as.matrix(sm)
   }

   if (is.null(w)){
      w <- rep(1,nrow(sm))
   }

   # wsm <- similarity_matrix_by_weights(sm, w)
   # return( mean( as.dist(wsm) ) )


   wlen <- length(w)
   wsum <- sum(w)


      sum_weighted_pairwise <- 0
      for (i in 1:wlen) {
         for (j in 1:wlen) {
            if (i > j) {
               sum_weighted_pairwise <- sum_weighted_pairwise + sm[i,j] * w[i] * w[j]
            }
            if (i == j) {
               sum_weighted_pairwise <- sum_weighted_pairwise + sm[i,j] * (w[i]-1) * w[j] / 2
            }

         }
      }

      number_weighted_pairwise <- wsum * (wsum-1) / 2
      mean_weighted_pairwise   <- sum_weighted_pairwise / number_weighted_pairwise

   return(mean_weighted_pairwise)
}

