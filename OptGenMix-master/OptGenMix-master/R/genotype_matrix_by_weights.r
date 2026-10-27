#' Transforms a genotype matrix by a vector of weights 
#' 
#' @param gt                - genotype matrix (samples by loci)
#' @param w                 - weights
#' @return gt matrix transformed according to weights 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export

genotype_matrix_by_weights <- function(gt, w) {

   n <- length(w)

   i_t <- NULL

   for (i in 1:n) {

      if (w[i] > 0) {

         iit <- rep(i, w[i])

         if (is.null(i_t)) {
            i_t <- iit
         } else {
            i_t <- c(i_t, iit)
         }

      } 
   } 

   gt_w <- gt[ i_t, ]

   return(gt_w)
}


