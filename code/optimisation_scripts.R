## Scripts from Jason Bragg
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


optimize_locations <- function(nsteps, sg, sd, gd ) {
  
  chain  <- mat.or.vec(nsteps,3)
  sols   <- list()
  
  csg    <- sg 
  cout   <- estimate_spatial_genetic_distances(sd, gd, csg)
  csum   <- summ_dists(cout,sg)   
  
  for (s in 1:nsteps) {
    
    psg  <- perform_switch(csg, i=csum$imax) 
    
    if (runif(1) > 0.5) {
      psg  <- perform_switch(psg) 
    }      
    
    pout <- estimate_spatial_genetic_distances(sd, gd, psg)
    psum <- summ_dists(pout,psg)   
    
    if (psum$imean < csum$imean) {
      csg  <- psg
      csum <- psum 
      chain[s,3] <- 1
    }
    
    chain[s,1] <- csum$imean
    chain[s,2] <- psum$imean
    sols[[s]]  <- csg
    
    cat(s, csum$imean,"\n")
  }
  
  return( list(chain=chain, sols=sols) )
}


optimize_locations_sa <- function(nsteps, sg, sd, gd, max_t ) {
  
  chain  <- mat.or.vec(nsteps,3)
  sols   <- list()
  
  csg    <- sg 
  cout   <- estimate_spatial_genetic_distances(sd, gd, csg)
  csum   <- summ_dists(cout,sg)   
  
  for (s in 1:nsteps) {
    
    psg  <- perform_switch(csg) 
    
    #      if (runif(1) > 0.5) {
    #         psg  <- perform_switch(psg, i=csum$imax) 
    #      }      
    
    
    pout <- estimate_spatial_genetic_distances(sd, gd, psg)
    psum <- summ_dists(pout,psg)   
    
    accept=FALSE
    
    if (psum$imean < csum$imean) {
      accept=TRUE
    } else {
      
      proportion_ahead = 1 - (s-1)/nsteps
      pa = proportion_ahead
      t <- pa * max_t
      
      delta <- psum$imean - csum$imean
      p_accept_worse <- exp( ( -delta ) / t  )
      
      if ( runif(1) < p_accept_worse ) {
        accept=TRUE
      }
    }
    
    
    if (accept) {
      csg  <- psg
      csum <- psum 
      chain[s,3] <- 1
    }
    
    chain[s,1] <- csum$imean
    chain[s,2] <- psum$imean
    sols[[s]]  <- csg
    
    cat(s, csum$imean,"\n")
  }
  
  return( list(chain=chain, sols=sols) )
}

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

summ_dists <- function(gdistlist, sg) {
  
  pmax      <- which.max((1/gdistlist$spatial_genetic[,1])*gdistlist$spatial_genetic[,2])
  imean     <- mean((1/gdistlist$spatial_genetic[,1])*gdistlist$spatial_genetic[,2],na.rm=TRUE)
  rnmax     <- sample(c(gdistlist$position_sample[pmax,1],gdistlist$position_sample[pmax,2]))[1]
  imax      <- which(names(sg) == rnmax)
  
  return(list(pmax=pmax, imean=imean, rnmax=rnmax, imax=imax))
}

