# R

### exploration

## attach resistance scores to inferred pedigree
rustScore <- meta$Score
names(rustScore) <- meta$accession

### 1 = resistant; 2 = susceptible
rustResRes <- meta$Score
names(rustResRes) <- meta$accession
rustResRes[meta$Score == 3 ] <- 1 

rustScore[rustScore == 2 ] <- 1
rustScore[rustScore == 3 ] <- 2

seedling_pedigree_rust_score <- cbind(seedling_pedigree, rustScore[ seedling_pedigree[,2]  ],rustResRes[ seedling_pedigree[,2]  ])

table(seedling_pedigree_rust_score[,7],seedling_pedigree_rust_score[,9])
table(seedling_pedigree_rust_score[,8],seedling_pedigree_rust_score[,9])

pkadisc <- pka
pkadisc[pka < 0.35] <- 0
pkadisc[pka >= 0.35] <- 1

# Later developed rust

alreadyRemoved <- "P2020-0371.S10"

nurseryInfections <- c(
"P2020-0371.S3",
"P2020-0373.S7",
"P2020-0246.S10"
)


inurseryInfections <- match (nurseryInfections, seedling_pedigree_rust_score[,2])
iunavail <- match(inurseryInfections, iavailable)
iavail <- iavailable[ -iunavail]

sprsa <- seedling_pedigree_rust_score[iavail,]
pkaa <- pk[sprsa[,1], sprsa[,1]]

spaa <- sharedParents[iavail, iavail]


#sanity check
cbind(rownames(pkaa), rownames(spaa), sprsa[,1:2])


locs_ABGMA <- 47
locs_Minyu <- 49

library(OptGenMix)
set.seed(98237)
Nt <- 96
max_wts <- rep(1, nrow(pkaa))


opt_pk_96  <- optimize_single_objective( sm=pkaa, N_t=Nt, weights_max=max_wts, measure="negative_matrix_weighted_mean", max_steps=50000, max_t=0.005, p_depends_delta=TRUE)


out_96 <- cbind(seedling_pedigree_rust_score[iavail,],opt_pk_96$weight[50000,])

table(out_96[,"MatLine"],out_96[,11])

only_96 <- out_96[which(out_96[,11]=="1"),]


table(only_96[seq(1,95,2),"MatLine"])

ABGMA_plants <- only_96[seq(1,95,2),]
ABGMA_plants_47 <- ABGMA_plants[-c(45),]
ABGMA_plants_47_names <- ABGMA_plants_47[,1]
pkAa <- pkaa[ABGMA_plants_47_names, ABGMA_plants_47_names]

x = 7
y = 7

nsites <- x*y
grid <- mat.or.vec(nsites, 2)
c <- 1
for (i in 1:x) {
   for (j in 1:y) {
         grid[c,] <- c(i, j)
         c <- c + 1
   }
}

grid <- grid[-c(1,49),]


location_names <- paste0("S",1:nrow(grid))
spatial_distance <- as.matrix(dist(grid))
colnames(spatial_distance) <- location_names
rownames(spatial_distance) <- location_names

location_samples <- sample(rownames(pkAa))
names(location_samples) <- location_names


source("/home/jgb/projects/hibb_planting/estimate_spatial_genetic_distances.r")
source("/home/jgb/projects/hibb_planting/perform_switch.r")
source("/home/jgb/projects/hibb_planting/optimize_locations.r")
source("/home/jgb/projects/hibb_planting/summ_dists.r")

nsteps <- 50000

set.seed(4589294)
optimized_location_chain_sa <- optimize_locations_sa(nsteps, location_samples,  spatial_distance, pkAa, 0.00001)

solution <- optimized_location_chain_sa$sols[[nsteps]]


save.image("/home/jgb/projects/rhodamnia/seedlingPedigreeOpt03Sept2024.Rda")

cbind(grid,names(solution),solution)

x = 7; y = 7
layout <- matrix(rep("",x*y),nrow=x)

for (i in 1:length(solution)) {
         l1 <- grid[i,1] 
         l2 <- grid[i,2] 
         layout[l1,l2] <- seedling_pedigree[which(seedling_pedigree[,1]==solution[i]),2]
}




write.csv(layout,"/home/jgb/projects/rhodamnia/ABGMA_planting_layout03Sept2024.csv")

iABGMArem <- which(only_96[,1] == "NSW1166484")
iminyu <- c(seq(1,95,2)+1,iABGMArem)

MINYU_plants <- only_96[ iminyu ,]
MINYU_plants_names <- MINYU_plants[,1]


### check
table(names(c(MINYU_plants_names, ABGMA_plants_47_names)))
length(table(names(c(MINYU_plants_names, ABGMA_plants_47_names))))
which(table(names(c(MINYU_plants_names, ABGMA_plants_47_names)))>1)


write.csv(cbind(rownames(MINYU_plants), MINYU_plants[,1]),"/home/jgb/projects/rhodamnia/MINYU_planting_list15Oct2024.csv", row.names=FALSE, col.names=FALSE, quote=FALSE)
