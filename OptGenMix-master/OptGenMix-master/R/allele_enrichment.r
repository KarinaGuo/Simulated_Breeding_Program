#' Calculates enrichment of desirable allele. Given genotypes (gt),
#' a vector (w) of weight values (equal length to the 
#' number of individuals, a vector (v) of values describing
#' the importance of each locus, estimates an index of enrichment of
#' preferred alleles. Assumes genotype fitness 2 > 1 > 0.  
#' 
#' @param gt  - genotype matrix for adaptive loci (fitness 2 > 1 > 0)
#' @param w   - vector of individual weights
#' @param loc - vector of locus adaptation values
#' @param rec - if true, fitness 2 > 1 = 0
#' @return An index of adaptedness, float value 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export

allele_enrichment <- function(gt, w=NULL, loc=NULL, rec=FALSE) {

   if (rec) {
      gt[ gt == 1 ] <- 0 
   }

   if (is.null(w)) {
      p  <- colSums(gt) / (nrow(gt) *2)
   } else {
      wm <- matrix(rep(w,ncol(gt)),ncol=ncol(gt),byrow=FALSE) 
      p  <- colSums(gt*wm) / (nrow(gt)*2) / mean(w)
   }

   if (is.null(loc)) {
      a  <- sum(p) / length(p)
   } else {
      a  <- sum(p*loc) / sum(loc)
   }

   return(a)
}

