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


optimize_multi_objective <- function( v1=NULL, v2=NULL, measure_1=NULL, measure_2=NULL, max_steps=10000, N_t=NULL, initial_weights=NULL, weights_max=NULL, weights_min=NULL, max_t=1, q=NULL, p_depends_delta=FALSE, disp=0, c1=1, c2=1, cboth=1, nda=FALSE, min_t=0, nd_samples=100) {

   cat( "\n\n" )
   cat( "  Multi-objective optimization commencing \n" )
   proceed=TRUE

   ### parse arguments
   # check for gt or sm
   if ( is.null(v1) ) {
      cat( "   Provide at least one matrix of values \n" )
      proceed=FALSE
   } else {
      N_g = nrow(v1) 
   }


   if ( is.null(v2) ) {
      cat( "   Both objective criteria will be applied to v1 \n" )
      v2=v1
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


   summary_1 <- NULL
   summary_2 <- NULL

   
   # generate a value of objective measure for initial
   summary_1 <- generate_measure_value(v1, measure=measure_1, w=initial_weights, q=q, disp=disp)
   summary_2 <- generate_measure_value(v2, measure=measure_2, w=initial_weights, q=q, disp=disp)

   # if objective measure initial returns NULL, problem
   if (is.null(summary)) {
      proceed=FALSE
      cat( "   Initial weights returned null value for measure \n" )
   }

   # if creating a non-dominated archive
   if (nda) {

      archive      <- list()
      nda_complete <- FALSE

   } else {

      nda_complete <- FALSE

   }



   # if proceed is still TRUE, we are ready to run
   if (!proceed) {
      cat( "   Something has gone wrong. Reconsider data or settings. \n" )
   } else {  

      ### run Simulated Annealing chain

      ### set up list to store chain
      weight  <- mat.or.vec(max_steps,N_g)
      value_1 <- mat.or.vec(max_steps,1)
      value_2 <- mat.or.vec(max_steps,1)
      accept  <- mat.or.vec(max_steps,1)
      chain   <- list(weight=weight, value_1=value_1, value_2=value_2, accept=accept)
      chain$value_1[1]   <- summary_1
      chain$value_2[1]   <- summary_2
      chain$weight[1,] <- initial_weights

      if (nda) {

         cat("   First sample added to non-dominated archive  \n")
         archive[[1]] <- list(value_1=summary_1, value_2=summary_2, weights=initial_weights)

      }


      # initialize some params to begin
      s <- 2
      weights <- initial_weights
      cat("   Commencing optimization \n")

      # MAIN CHAIN LOOP
      while ( s <= max_steps & !nda_complete ) {

         proposed_weights   <- propose_new_weights(weights, w_max=weights_max, w_min=weights_min)

         proposal_summary_1 <- generate_measure_value(v1, measure=measure_1, w=proposed_weights, q=q, disp=disp)
         proposal_summary_2 <- generate_measure_value(v2, measure=measure_2, w=proposed_weights, q=q, disp=disp)

         temp               <- temp_scheduler(s, max_steps, max_t=max_t)

         # hold at moderate temp to allow
         # exploration of nd solutions
         if (nda) {
            if (temp < min_t) {
               temp <- min_t
            }
         }

         accept_proposal    <- multi_accept_reject(summary_1=summary_1, summary_2=summary_2, 
                                                   proposal_summary_1=proposal_summary_1, proposal_summary_2=proposal_summary_2, 
                                                   temp, p_depends_delta=p_depends_delta, c1=c1, c2=c2, cboth=cboth)


         if (accept_proposal) {
            weights <- proposed_weights
            summary_1 <- proposal_summary_1
            summary_2 <- proposal_summary_2
            
            if (nda) {
               archive <- reconcile_sample_nondominated_archive(summary_1, summary_2, weights, archive)
               # cat("   Sample ", length(archive), " of ", nd_samples, " added to non-diminated archive  \n")
               if ( length(archive) >= nd_samples ) { nda_complete <- TRUE }
            }

         } 

         chain$accept[s]    <- accept_proposal
         chain$value_1[s]   <- summary_1
         chain$value_2[s]   <- summary_2
         chain$weight[s,]   <- weights

         s <- s + 1

         pfreq <- floor(max_steps / 20)
         if (s %% pfreq == 0) {cat("   Step:", s, "\n")}

      } # end MAIN CHAIN LOOP

      if (nda) {

         value_1 <- unlist(do.call("cbind", archive)[1,])
         value_2 <- unlist(do.call("cbind", archive)[2,])
         archive_values <- cbind(value_1, value_2)

         archive_weights          <- matrix( unlist(do.call("cbind", archive)[3,]),nrow=length(archive),byrow=TRUE)
         archive_matrices         <- list(archive_values=archive_values, archive_weights=archive_weights)
         chain[["archive"]]       <- archive_matrices
      }

      return(chain)

   } # end if proceed

}

