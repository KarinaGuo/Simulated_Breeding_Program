perform_switch <- function(sg, i=NULL, j=NULL) {

   ls <- sg
   len <- length(sg)

   if (is.null(i)) {
     i = sample(1:len)[1]
   }

   if (is.null(j)) {
     j = sample(1:len)[1]
   }

   ls[i] <- sg[j]
   ls[j] <- sg[i]

   return(ls)
}

