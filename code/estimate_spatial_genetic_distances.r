estimate_spatial_genetic_distances <- function(sd, gd, sg) {

nloc <- nrow(sd)

out_names <- mat.or.vec(nloc*(nloc-1)/2,4)
out_nums  <- mat.or.vec(nloc*(nloc-1)/2,2)

line <- 1
for (i in 1:nloc) {

   for (j in 1:nloc) {

      if (i < j) {
        is <- sg[i]
        js <- sg[j]

        spatial <- sd[ i, j ]

        if (is.na(is) | is.na(js)) {
           genetic <- NA 
        } else {
           genetic <- gd[ is, js ]
        }
        out_names[line,] <- c(rownames(sd)[i], rownames(sd)[j], is, js)
        out_nums[line,] <- c(spatial, genetic)
        line <- line + 1
      }
   }
}

return(list(position_sample=out_names, spatial_genetic=out_nums))

}

