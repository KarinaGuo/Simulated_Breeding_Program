#' Calculates a mean value for a nominated vector (vs), given
#' a set of weights (w). Returns the weighted mean of the absolute 
#' values of differences between vector values and a 'displacement' 
#' This could be used, e.g., to minimize the mean difference between  
#' temperature of origin for each sample, and temperate of a site.  
#' 
#' @param vs   - genotype matrix (individuals, loci)
#' @param w    - vector of weights
#' @param disp - opt value -- displacement from zero
#' @return difference 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export


sum_of_squared_difference <- function(vs, w, disp=0) {

   ss <- sum( ((vs-disp)^2) *w) 
   return(ss)

}
