#' Calculates a mean value for a nominated vector (vs), given
#' a set of weights (w). Returns the absolute value of the 
#' difference between this weighted mean, and a nominated value
#' (diff). This could be used to minimize the use of individuals 
#' low levels of heterozygosity, or to place arbitrary 
#' constraints on genotype composition 
#' 
#' @param vs   - genotype matrix (individuals, loci)
#' @param w    - vector of weights
#' @param disp - opt value -- displacement from zero
#' @return difference 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export


weighted_mean_of_vector <- function(vs, w, disp=0) {

   ss <- abs((sum(vs*w) / sum(w)) - disp)
   return(ss)

}
