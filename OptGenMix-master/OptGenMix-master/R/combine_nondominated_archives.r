#' Combines two non-dominated archives, output from 
#' multiobjective optimization function
#' 
#' @param a1                - a non-dominated archive
#' @param a2                - a non-dominated archive
#' @return a list of output 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export

combine_nondominated_archives <- function(a1, a2) {

   anew   <- a1

   anew$archive_values <- rbind(anew$archive_values, a2$archive_values)
   anew$archive_weights <- rbind(anew$archive_weights, a2$archive_weights)

   anew_len <- nrow(anew$archive_values)

   i_dom <- NULL

   for (i in 1:anew_len) {

      iv1 <- anew$archive_values[i,1]
      iv2 <- anew$archive_values[i,2]


      for (j in 1:anew_len) {

         iv1_smaller <- FALSE
         iv2_smaller <- FALSE 

         if ( i != j ) { 
            jv1 <- anew$archive_values[j,1]
            jv2 <- anew$archive_values[j,2]

            if (iv1 < jv1) {iv1_smaller <- TRUE}
            if (iv2 < jv2) {iv2_smaller <- TRUE}

            if ( iv1_smaller & iv2_smaller ) { 
                  i_dom <- c(i_dom, i)
                  i_dom <- unique(i_dom) 
            }
                     
         }

      }

   }

   if (length(i_dom) > 0) {    
      anew$archive_values  <- anew$archive_values[-i_dom,]
      anew$archive_weights <- anew$archive_weights[-i_dom,]
   } 
   return(anew)

}


