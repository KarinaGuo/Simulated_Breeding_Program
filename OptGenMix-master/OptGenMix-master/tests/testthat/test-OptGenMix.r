library(testthat)
library(OptGenMix)


################################
# Shannon tests based on values
# in Sherwin et al. TREE article

test_that("shannon_TREE_examples_fixed",{
  gt0_1 <- matrix(c(0,0),ncol=1) 
  s0_1 <- shannon_diversity(gt0_1,q=0)

  gt0_2 <- matrix(c(2,2),ncol=1) 
  s0_2 <- shannon_diversity(gt0_1,q=0)

  s1_1 <- shannon_diversity(gt0_1,q=1)
  s2_1 <- shannon_diversity(gt0_1,q=2)

  expect_equal( s0_1, s0_2)
  expect_equal( s0_1, 1)
  expect_equal( s1_1, 1)
  expect_equal( s2_1, 1)
})


test_that("shannon_TREE_examples_quarter",{
  gt0 <- matrix(c(1,0),ncol=1) 
  s0  <- shannon_diversity(gt0,q=0)
  s1  <- shannon_diversity(gt0,q=1)
  s2  <- shannon_diversity(gt0,q=2)

  expect_equal( s0, 2)
  expect_equal( s1, exp(0.56), tolerance=0.02)
  expect_equal( s2, 1/(1-0.38),tolerance=0.02)

})

test_that("shannon_TREE_examples_half",{
  gt0 <- matrix(c(1,1),ncol=1) 
  s0  <- shannon_diversity(gt0,q=0)
  s1  <- shannon_diversity(gt0,q=1)
  s2  <- shannon_diversity(gt0,q=2)

  expect_equal( s0, 2)
  expect_equal( s1, exp(0.69), tolerance=0.02)
  expect_equal( s2, 1/(1-0.5),tolerance=0.02)

})


################################
# Nei test -- expected Heterozygosity

test_that("Nei_example",{
  gt  <- rbind(c(2,2,0,0), c(1,2,2,1), c(2,2,1,0), c(1,2,2,0), c(2,1,0,0), c(2,0,0,1)) 
  nei <-  nei_diversity(gt)

  p <- colSums(gt) / nrow(gt) / 2  
  Hexp <- 2*p*(1-p)

  expect_equal( nei, mean(Hexp))
})



################################
# Allele enrichment tests

test_that("allele_enrichment",{

  gt  <- rbind(c(2,2,0,0), c(1,2,2,1), c(2,2,1,0), c(1,2,2,0), c(2,1,0,0), c(2,0,0,1)) 

  adp <-  allele_enrichment(gt)
  adpv <-  sum(gt) / (nrow(gt)*ncol(gt)*2)

  expect_equal( adp, adpv)

  adph <-  allele_enrichment(gt, rec=TRUE)
  gth  <- gt
  gth[gt == 1]  <- 0 
  adphv <-  sum(gth) / (nrow(gth)*ncol(gth)*2)
  expect_equal( adph, adphv)

  lvec <- c(0,1,1,0)

  adp_lvec <-  allele_enrichment(gt, loc=lvec)
  adpv_lvec <-  sum(gt[,2:3]) / (nrow(gt[,2:3])*ncol(gt[,2:3])*2)
  expect_equal( adp_lvec, adpv_lvec)

})


# test weighted matrix mean

test_that("weighted_mean_matrix",{

   sm <- rbind( c(1,0,0), c(0,1,0), c(0,1,1))
   w  <- c(1,2,2)

   wms <- weighted_mean_of_pairwise_matrix(sm, w)
   expect_equal( wms, 0.6)


})



# test weighted mean vector

test_that("weighted_mean_vector",{

   vc <- c(0,0,1)
   w  <- c(1,2,2)


   wms <- weighted_mean_of_vector(vc, w)
   expect_equal( wms, 0.4)

})


test_that("sum_of_squared_difference",{

   v1 <- c(1,2,1,2)
   v2 <- c(1,2,1,6)
   w1 <- c(1,1,1,1)
   w2 <- c(0,1,0,1)
   w3 <- c(1,1,1,3)
   d1 <- 1.5
   d2 <- 2

   ssd1 <- sum_of_squared_difference(v=v1,w=w1,disp=d1)
   expect_equal( ssd1, 1)

   ssd2 <- sum_of_squared_difference(v=v1,w=w2,disp=d2)
   expect_equal( ssd2, 0)

   ssd3 <- sum_of_squared_difference(v=v2,w=w3,disp=d2)
   expect_equal( ssd3, 50)

})


test_that("temp_1",{
  s     <- 1
  smax  <- 1000
  max_t <- 100
  t1 <- temp_scheduler(s,smax, max_t=max_t)

  expect_equal(t1, (max_t - (max_t/smax)) )
})


test_that("temp_smax",{
  s     <- 1000
  smax  <- 1000
  max_t <- 100
  t1000 <- temp_scheduler(s,smax, max_t=max_t)

  expect_equal(t1000, 0 )
})



test_that("intial_weight_random",{
  N_g <- 100
  N_t <- 1000

  p <- propose_initial_weights(N_g, N_t)
  expect_equal(mean(p), 10, tolerance=2 )
})


test_that("intial_weight_constraint",{
  N_g <- 100
  N_t <- 1000

  w_max <- rep(N_t,N_g)
  w_max[1:50] <- 0

  p <- propose_initial_weights(N_g, N_t, w_max=w_max)
  expect_equal(mean(p[51:100]), 20, tolerance=3 )
})


test_that("propose_new_weights",{
   N_g <- 10
   N_t <- 100
   wold <- rep(N_t / N_g, N_g)
   n <- 10000
   out <- mat.or.vec(n, N_g)

   add <- mat.or.vec(n, N_g)
   rem <- mat.or.vec(n, N_g)

   add[ out==11] <- 1
   rem[ out==9] <- 1
   for (i in 1:n) {
      out[i,] <- propose_new_weights(wold)
   }
   add[ out==11] <- 1
   rem[ out==9] <- 1

   nadd <- colSums(add)
   nrem <- colSums(rem)

   expect_equal( max(nadd), n/N_g, tolerance=0.1, scale=n/N_g)
   expect_equal( min(nadd), n/N_g, tolerance=0.1, scale=n/N_g)
   expect_equal( max(nrem), n/N_g, tolerance=0.1, scale=n/N_g)
   expect_equal( min(nrem), n/N_g, tolerance=0.1, scale=n/N_g)

})


test_that("generate_measure_value",{

  gt  <- rbind(c(2,2,0,0), c(1,2,2,1), c(2,2,1,0), c(1,2,2,0), c(2,1,0,0), c(2,0,0,1)) 
  nei <-  nei_diversity(gt)
  gmv <- generate_measure_value(v=gt, measure="nei")
  expect_equal( nei, gmv)

  w <- c(0,2,2,2,0,0)
  wnei <-  nei_diversity(gt,w=w)
  wgmv <- generate_measure_value(v=gt, w=w, measure="nei")
  expect_equal( wnei, wgmv)

  w <- c(0,2,2,2,0,0)
  wshan <-  shannon_diversity(gt,w=w,q=2)
  wwshan <- generate_measure_value(v=gt, w=w, measure="shannon", q=2)
  
})


# archive <- list()
# archive[[1]] <- list(value_1=1, value_2=1, weights=c(1,1,1,1,1))
# archive <- reconcile_sample_nondominated_archive(1.5,0.5,c(0,1,1,3),archive)


#disp <- 0.5
#w <- c(1,1,1,1,1,1)


#vs <- c(0,1,0,1,0,1)
#abs((sum(vs*w) / sum(w)) - disp)
#sum(abs(vs-disp)*w) / sum(w)
#weighted_mean_of_absolute_difference(vs, w, disp)
#weighted_mean_of_vector(vs,w,disp)

#vs <- c(0.25,0.75,0.25,0.75,0.25,0.75)
#abs((sum(vs*w) / sum(w)) - disp)
#sum(abs(vs-disp)*w) / sum(w)

#weighted_mean_of_absolute_difference(vs, w, disp)
#weighted_mean_of_vector(vs,w, disp)
