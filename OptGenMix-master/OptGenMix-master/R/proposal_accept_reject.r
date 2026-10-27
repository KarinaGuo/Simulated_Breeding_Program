#' Determines acceptance or rejection of a proposal
#' 
#' @param summary          - summary of current
#' @param proposal_summary - summary of new proposal
#' @param t                - current temperature
#' @param p_depends_delta  - if TRUE, accept depends on difference between new and old 
#' @param c                - constant multiplier for numerator (default, c=1)
#' @return Boolean, true if proposal accepted 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export

proposal_accept_reject <- function(summary=summary, proposal_summary=proposal_summary, t, p_depends_delta=FALSE, c=1) {

   accept=FALSE

   # if better, accept
   if ( proposal_summary > summary ) {

      accept=TRUE

   } else {

      delta <- summary - proposal_summary 
      if (p_depends_delta) {
         numerator <- delta * c
      } else {
         numerator <- c
      }

      p_accept_worse <- exp( ( -numerator ) / t  )

      #cat(p_accept_worse, "\n")
      if ( runif(1) < p_accept_worse ) {
         accept=TRUE
      }
   }

   return(accept)
}


