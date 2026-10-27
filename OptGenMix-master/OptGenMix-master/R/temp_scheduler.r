#' Determines the temperature for annealing
#' 
#' @param s         - current time step
#' @param max_steps - maximum step
#' @return temperature, floating point value 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export

temp_scheduler <- function(s, max_steps, max_t) {

   proportion_ahead = 1 - s/max_steps
   pa = proportion_ahead

   t <- pa * max_t

   return(t)
}
