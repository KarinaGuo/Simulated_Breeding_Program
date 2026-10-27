#' Performs simulated annealing to find mixture of individuals 
#' that maximize a function (e.g. a measure of genetic diversity) 
#' 
#' @param gt                - a genotype matrix (samples, loci)
#' @param initial_weights   - a vector of initial weights
#' @param weights_max       - maximum weights for each individual
#' @param weights_min       - minimum weights for each individual
#' @param max_t             - maximum temperature 
#' @return a list of output 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export


optimize_single_objective <- function( gt=NULL, sm=NULL, measure=NULL, max_steps=10000, N_t=NULL, initial_weights=NULL, weights_max=NULL, weights_min=NULL, max_t=1, q=NULL, p_depends_delta=FALSE, disp=0) {

   proceed=TRUE

   ### parse arguments
   # check for gt or sm
   if ( is.null(gt) & is.null(sm)) {
      cat( "   Provide either initial_weights or a number of target individuals (N_t) \n" )
      proceed=FALSE
   }

   if ( !is.null(gt) & !is.null(sm)) {
         cat( "   Provide either a genotype or dimilarity matrix \n" )
         proceed=FALSE
   }

   # check for availability of N_g and store value
   N_g=NULL

   if (!is.null(gt)) {
      N_g = nrow(gt) 
      v = gt
   }

   if (!is.null(sm)) {
      N_g = nrow(sm)
      v = sm 
   }

   # check for N_t or initial_weights
   if ( is.null(N_t) & is.null(initial_weights)) {
      cat( "   Provide either initial_weights or a number of target individuals (N_t) \n" )
      proceed=FALSE
   } 

   if ( !is.null(N_t) & !is.null(initial_weights)) {
      if (N_t != sum(initial_weights)) {
         cat( "   Conflict in number of target individuals and initial weights provided \n" )
         proceed=FALSE
      }
   }



   # if N_t but no initial weights, generate initial weights
   if ( !is.null(N_t) & is.null(initial_weights)) {
      cat( "   Generating initial weights  \n" )
      initial_weights <- propose_initial_weights(N_g, N_t)
   }



   summary <- NULL

   
   # generate a value of objective measure for initial
   summary <- generate_measure_value(v, measure=measure, w=initial_weights, q=q, disp=disp)

   # if objective measure initial returns NULL, problem
   if (is.null(summary)) {
      proceed=FALSE
      cat( "   Initial weights returned null value for measure \n" )
   }


   # if proceed is still TRUE, we are ready to run
   if (!proceed) {
      cat( "   Something has gone wrong. Reconsider data or settings. \n" )
   } else {  

      ### run Simulated Annealing chain

      ### set up list to store chain
      weight  <- mat.or.vec(max_steps,N_g)
      value   <- mat.or.vec(max_steps,1)
      accept  <- mat.or.vec(max_steps,1)
      chain   <- list(weight=weight, value=value, accept=accept)
      chain$value[1]   <- summary
      chain$weight[1,] <- initial_weights

      # initialize some params to begin
      s <- 2
      weights <- initial_weights

      # MAIN CHAIN LOOP
      while ( s <= max_steps ) {

         proposed_weights <- propose_new_weights(weights, w_max=weights_max, w_min=weights_min)

         proposal_summary <- generate_measure_value(v, measure=measure, w=proposed_weights, q=q, disp=disp)

         temp             <- temp_scheduler(s, max_steps, max_t=max_t)

         accept_proposal  <- proposal_accept_reject(summary=summary, proposal_summary=proposal_summary, temp, p_depends_delta=p_depends_delta)
   
         if (accept_proposal) {
            weights <- proposed_weights
            summary <- proposal_summary
         } 

         chain$accept[s]  <- accept_proposal
         chain$value[s]   <- summary
         chain$weight[s,] <- weights

         s <- s + 1

         pfreq <- floor(max_steps / 20)
         if (s %% pfreq == 0) {cat("   up to step", s, "\n")}

      } # end MAIN CHAIN LOOP

      return(chain)

   } # end if proceed

}

