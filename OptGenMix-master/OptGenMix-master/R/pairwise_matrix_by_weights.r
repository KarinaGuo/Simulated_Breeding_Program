#' Transforms a pairwise matrix by a vector of weights 
#' 
#' @param sm                - pairwise sample matrix (samples by samples)
#' @param w                 - weights
#' @return sm matrix transformed according to weights 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export

pairwise_matrix_by_weights <- function(sm, w) {

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

   sm_w <- sm[ i_t, i_t]

   return(sm_w)
}

