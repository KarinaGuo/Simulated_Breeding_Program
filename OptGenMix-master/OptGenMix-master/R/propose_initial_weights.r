#' Proposes initial weights given number of genotypes (length of weights) 
#' and number of individuals in target population (sum of weights).   
#' 
#' @param N_g    - number of individuals with genotypes
#' @param N_t    - number of individuals in target population
#' @param w_max  - a vector of maxima for each individual
#' @return An initial, random, set of weights
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export


propose_initial_weights <- function(N_g, N_t, w_max=NULL, verbose=FALSE) {

   if (verbose) {cat("   Generating initial weights at random \n")}

   # get an initial set of weights by sampling with replacement
   vector_initial    <- sample(N_g,size=N_t,replace=TRUE)

   # simpler to use table (but in case re-implemented)
   collapse_weights <- function(vi, N_g) {

      wvec <- rep(0, N_g)
      for (i in 1:N_g) {
         i_count <- length( which(vi==i) )  

         if (i_count > 0) {
            wvec[i] <- i_count
         } 
      }
      return(wvec)
   }

   weights_initial <- collapse_weights(vi=vector_initial, N_g=N_g)

   
   if ( is.null(w_max) ) {
 
      if (verbose) {cat("   No constraints specified -- returning values \n")}
      w_out <- weights_initial

   } else {

      if ( sum(w_max) < N_t ) {
         cat("  Sum of constraint values is smaller than total target population \n")
         cat("  Returning initial random weights... \n")
         cat("  Reconsider maximum weight vector... \n")
         return(weights_initial)
      }

      # now, if max weights are supplied, we will come up with an  
      # initial solution conforming to those values   
      if ( !(any(weights_initial > w_max)) ) {

         w_out <- weights_initial
         if (verbose) { cat("   Random weights satisfy constraints -- returning values \n") }

      } else {

         max_iter <- 100
         iters    <- 0
         w <- weights_initial

         if (verbose) {
            cat("   Initial random weights inconsistent with constraints \n")
            cat("   Generating new values  \n")
         }

         while( any(w > w_max) & iters < max_iter ) {

            ind_max_ge <- which(w >= w_max)
            sum_max_ge <- sum(w[ind_max_ge] - w_max[ind_max_ge])

            pinc <- rep(1,N_g)
            pinc[ind_max_ge] <- 0

            # choose some genotypes to add
            i_add  <- sample(N_g,size=sum_max_ge,replace=TRUE,pinc)
            add    <- collapse_weights(i_add,N_g)

            # find genotypes above maximum
            rem    <- rep(0, N_g)
            rem[ind_max_ge] <- w[ind_max_ge] - w_max[ind_max_ge]

            # remove genotypes
            w <- w - rem
         
            # add genotypes
            w <- w + add         

            # don't let it keep going if no solution
            iters <- iters + 1

            if (any(w > w_max) & iters >= max_iter) {
               if (verbose) {
                  cat("   Constraints not satisfied after", max_iter, " iterations...\n")
                  cat("   Consider a manual formulation of constraint values. \n")
               }
            }         
 
            if ( !(any(w > w_max)) ) {
               if (verbose) {
                  cat("  Constraints satisfied after", iters, " iterations of searching.\n")
               }
            }

         } # while w > w_max

         w_out <- w

      } # if any()

   } # if else is.null(w_max)

   return(w_out)
}

