######################################
# Helper Functions
######################################
#take in a vector of breeding values (i.e. the genetic component of trait w/ no GxE)
#sorted by population, the number of individuals in each population, and ploidy.
#return an the "RB" (Relethford-Blangero) Qst estimate using the Weaver 2017 notation.
Qst.fromtraitvals.rb <- function(bv, n, ploidy){
  if(length(n)==1){n <- rep(n, length(bv)/n)}
  popvec <- rep(1:length(n), times = n)
  grp.var <- aggregate(bv, by = list(popvec), function(x){var(x)})[,2]
  win.var <- sum(grp.var*n)/sum(n)
  btw.var <- var(aggregate(bv, by = list(popvec), function(x){mean(x)})[,2])*(length(n) - 1)/length(n)
  as.numeric(btw.var / (btw.var + 2*win.var))
}

popvec <- rep(1:length(10), times = 10)
bv <- rnorm(10)
grp.var <- aggregate(bv, by = list(popvec), function(x){var(x)})[,2]
sum(var(bv))

#return a "PBS" (Prout-Barker-Spitze) Qst using Weaver 2017 notation with identical input
#this is Qst estimated using the unbiased variance estimates
Qst.fromtraitvals.pbs <- function(bv, n, ploidy){
  if(length(n)==1){n <- rep(n, length(bv)/n)}
  popvec <- rep(1:length(n), times = n)
  grp.var <- aggregate(bv, by = list(popvec), function(x){var(x)})[,2]
  win.var <- sum(grp.var*n)/sum(n)
  btw.var <- var(aggregate(bv, by = list(popvec), function(x){mean(x)})[,2])
  as.numeric(btw.var / (btw.var + 2*win.var))
}

#generate one vector of random breeding values
#nloci must be no greater than the number of loci in the genotype matrix
phens.nloci <- function(genmat, n.loci, es.dist, ploidy=2, ...){
  loci.pick <- sample(1:nrow(genmat))[1:n.loci]
  if(is.null(es.dist)){ #alpha model with alpha = -1
    p <- mean(genmat[loci.pick,])/ploidy
    dis.std <- sqrt(1/(2*p*(1-p)))
    betas <- rnorm(n.loci, 0, dis.std)
  }else{
    betas <- es.dist(n.loci, ...)
  }
  if(n.loci == 1){
    return(genmat[loci.pick,] * betas)
  }
  as.numeric(t(genmat[loci.pick,])%*%betas)
}

#get an empirical distribution of Qsts given a genotype matrix,
#a number of causal loci (chosen from among the genotypes, must be no greater
#than number of rows in genotype matrix), a number of traits to simulate,
#vector of number of individuals per population, distribution from which to draw
#effect sizes, a ploidy, and any extra parameters for the effect size distribution
Qst.dists <- function(genmat, n.loci, n.trait, n, es.dist, ploidy= 2, ...){
  Qsts.rb <- numeric(n.trait)
  Qsts.pbs <- numeric(n.trait)
  for(i in 1:n.trait){
    phens <- phens.nloci(genmat, n.loci, es.dist, ...)
    Qsts.rb[i] <- Qst.fromtraitvals.rb(phens, n, ploidy)
    Qsts.pbs[i] <- Qst.fromtraitvals.pbs(phens, n, ploidy)
  }
  cbind(Qsts.rb, Qsts.pbs)
}




######################################
# Simulation starts here
######################################
setwd("/Users/BidYD/Desktop/F/fst")
#install.packages("extraDistr")
library(extraDistr)

FST = 0.1
# dist = rnorm
# dist.name = "gaussian"

n.trait = 10000
inds = 2000

for (model in c("split_star", "split_balanced", "split_caterpillar", "migration_island", "migration_circular")){
  print(model)
  for (demes in c(2,4,8,16)){
    if (any(file.exists(paste0("input/fst", FST, "/geno_", demes, "D_", model, ".txt")))){
      print(paste0(demes, " demes"))
      gens <- data.matrix(read.table(paste0("input/fst", FST, "/geno_", demes, "D_", model, ".txt"), sep = ""))
      Qstmat <- matrix(ncol = 8, nrow = n.trait)
      n <- rep(inds, demes)
      for(j in c(1,3,5,7)){
        print(j)
        Qstmat[,j:(j+1)] <- Qst.dists(gens, n.loci = 10^((j-1)/2), n.trait = n.trait, n = n, es.dist = dist)
      }
      saveRDS(Qstmat, file.path(paste0("output/", dist.name, "/", demes, "D_", model, "_", FST ,"_qst_mtx.rds")))
    }
  }
}


