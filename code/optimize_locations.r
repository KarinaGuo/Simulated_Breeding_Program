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

