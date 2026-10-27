summ_dists <- function(gdistlist, sg) {

  pmax      <- which.max((1/gdistlist$spatial_genetic[,1])*gdistlist$spatial_genetic[,2])
  imean     <- mean((1/gdistlist$spatial_genetic[,1])*gdistlist$spatial_genetic[,2],na.rm=TRUE)
  rnmax     <- sample(c(gdistlist$position_sample[pmax,1],gdistlist$position_sample[pmax,2]))[1]
  imax      <- which(names(sg) == rnmax)

  return(list(pmax=pmax, imean=imean, rnmax=rnmax, imax=imax))
}

