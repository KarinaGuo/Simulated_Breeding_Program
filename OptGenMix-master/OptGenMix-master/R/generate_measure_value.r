#' Wraps functions that are used to calculate summaries 
#' 
#' @param v                 - gt or sm
#' @param measure           - "nei" or "shannon" or "mean_similarity"
#' @param w                 - weights
#' @param q                 - q value for shannon diversity 
#' @return summary value (if similarity, summary value * -1) 
#' @author Jason Bragg (jasongbragg@gmail.com)
#' @export

generate_measure_value <- function(v, measure=NULL, w=NULL, q=NULL, loc=NULL, disp=0) {

     if (is.null(measure)) {
        cat("   A measure is required \n")
        return(NULL)
     }

      if (measure=="nei") { 
            return( nei_diversity(gt=v, w=w) ) 
      }

      if (measure=="shannon") {
 
            if (is.null(q)) { 
               return( shannon_diversity(gt=v, w=w) ) 
            } else {

               if ( q==0 | q==1 | q==2) { 
                  return( shannon_diversity(gt=v, w=w, q=q ) ) 
               } else {
                  return(NULL)
               } 
            }
         
      }

      if (measure=="allele_enrichment") {
            return( 2*allele_enrichment(gt=v, w=w, loc=loc, rec=FALSE) )
      }

      if (measure=="negative_vector_weighted_mean") {
            return( -1 * weighted_mean_of_vector(vs=v, w=w, disp=disp) )
      }

      if (measure=="sum_squared_difference") {
            return( sum_of_squared_difference(vs=v, w=w, disp=disp) )
      }


      if (measure=="negative_sum_squared_difference") {
            return( -1 * sum_of_squared_difference(vs=v, w=w, disp=disp) )
      }


      if (measure=="vector_weighted_mean") {
            return( weighted_mean_of_vector(vs=v, w=w, disp=disp) )
      }

      if (measure=="negative_vector_diff_weighted_mean") {
            return( -1 * weighted_mean_of_absolute_difference(vs=v, w=w, disp=disp) )
      }

      if (measure=="vector_diff_weighted_mean") {
            return( weighted_mean_of_absolute_difference(vs=v, w=w, disp=disp) )
      }

      if (measure=="negative_matrix_weighted_mean") {
            return( -1*weighted_mean_of_pairwise_matrix(sm=v, w=w) )
      }

      if (measure=="matrix_weighted_mean") {
            return( weighted_mean_of_pairwise_matrix(sm=v, w=w) )
      }


}

