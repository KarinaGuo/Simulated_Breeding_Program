#' Calculates Shannon diversity, given a set of genotypes (gt)
#' and a vector of weight values (equal length to the 
#' number of genotyped individuals) 
#' 
#' @param gt  - genotype matrix (individuals, loci)
#' @param w   - vector of weights
#' @param q   - Shannon level (0, 1, or 2)
#' @return Shannon diversity, a float value 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export

shannon_diversity <- function(gt, w=NULL, q=1) {

   # caclulate allele frequencies
   if (is.null(w)) {
      p  <- colSums(gt) / (nrow(gt) *2)
   } else {
      wm <- matrix(rep(w,ncol(gt)),ncol=ncol(gt),byrow=FALSE) 
      p  <- colSums(gt*wm) / (nrow(gt)*2) / mean(w)
   }

   # q=0, allele count 
   if (q == 0) {

      a   <- p
      a[ which(p == 1) ] <- 1
      a[ which(p == 0) ] <- 1
      a[ which(p > 0 & p < 1) ] <- 2

      D_0_mean <- mean(a)

      return(D_0_mean)
   }

   # q=1, Shannon diversity 
   if (q == 1) {

      H_1 <- rep(0,ncol(gt))

      for (i in 1:ncol(gt)) {

         p_i <- p[i]
         q_i <- 1 - p_i 

         if (p_i == 0 | p_i == 1) {
            # 0 times log 0 has limit of 0
            # but R will return NaN
            H_1[i] <- - log(1) 

         } else {
            H_1[i] <- -1 * ( p_i * log(p_i) + q_i * log(q_i) ) 
         }

      }

      D_1 <- exp(H_1)

      D_1_mean <- mean(D_1) 

      return(D_1_mean)

   }

   # q=2 
   if (q == 2) {

      H_2 <- rep(0,ncol(gt))

      for (i in 1:ncol(gt)) {

         p_i <- p[i]
         q_i <- 1 - p_i 

         H_2[i] <- 1 - ( p_i*p_i + q_i*q_i) 

      }

      D_2 <- 1 / ( 1 - H_2 )

      D_2_mean <- mean(D_2) 

      return(D_2_mean)
   }


}

