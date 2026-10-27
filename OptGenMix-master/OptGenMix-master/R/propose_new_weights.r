#' Proposes new weights in the neighborhood of current 
#' state. Uses current state, constraints (maxima). 
#' 
#' @param w      - a vector of weights
#' @param w_max  - a vector of maxima for each individual
#' @param w_min  - a vector of minima for each individual
#' @return An updated set of weights 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export


propose_new_weights <- function(w, w_max=NULL, w_min=NULL) {

   # remove a single weight
   rem_weight_prob <- w

   # if min weights, apply
   if (!is.null(w_min)) {

      if ( any( w <= w_min ) ) {

         ind_le_min <- which(w <= w_min)
         rem_weight_prob[ ind_le_min ] <- 0 

      }

   }

   i_remove <- sample(1:length(w),1,FALSE, rem_weight_prob)
   w[i_remove] <- w[i_remove] - 1

   # initial weights all equal
   add_weight_prob <- rep(1, length(w))

   # don't add to individual just removed
   add_weight_prob[i_remove] <- 0 

   # if max weights, apply
   if (!is.null(w_max)) {

      if ( any( w >= w_max ) ) {

         ind_ge_max <- which(w >= w_max)
         add_weight_prob[ ind_ge_max ] <- 0 

      }

   }

   i_add <- sample(1:length(w),1,FALSE, add_weight_prob)
   w[i_add]   <- w[i_add] + 1

   return(w)
}

