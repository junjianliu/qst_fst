######################################
# Helper Functions
######################################
#take in matrix of genotypes (one row per locus), vector of sample sizes per pop 
#(or a scalar if sizes are identical), and ploidy (default is 2)
#return a vector with sample allele frequencies in the total population
sample.af.t <- function(genmat, ploidy=2){
  rowMeans(genmat)/ploidy
}

#return a matrix with sample allele frequencies in each subpopulation
sample.af.w <- function(genmat, n, ploidy=2){
  if(length(n)==1){n <- rep(n, ncol(genmat)/n)}
  popvec <- rep(1:length(n), times = n)
  af <- sample.af.t(genmat[,popvec==1], ploidy)
  for(j in 2:length(n)){
    af <- cbind(af, sample.af.t(genmat[,popvec==j], ploidy))
  }
  af
}

#return single-locus FST estimates via Nei's approach (Nei, 1973)
locus.Fsts.Nei <- function(genmat, n, ploidy=2){
  if(length(n)==1){n <- rep(n, ncol(genmat)/n)}
  af.w <- sample.af.w(genmat, n, ploidy)
  af.t <- sample.af.t(genmat, ploidy)
  ht <- af.t*(1-af.t)
  hw <- ((af.w*(1-af.w)) %*% n) / sum(n)
  (ht - hw)/ht
}

#return the single-locus FST estimates via the Weir-Cockerham approach 
#(with equal sample sizes, we followed eq. 7 in Weir and Cockerham, 1984)
locus.Fsts.WC <- function(genmat, n, ploidy=2){
  if(length(n)==1){n <- rep(n, ncol(genmat)/n)}
  n.samp <- sum(n)/length(n)
  r <- length(n)
  p.bar <- sample.af.t(genmat, 2)
  s.sqr <- apply(sample.af.w(genmat, n, 2), 1, function(x) var(x))
  h.bar <- apply(genmat, 1, function(x) sum(x==1)/sum(n))
  hb <- s.sqr - (p.bar*(1-p.bar) - (r-1)/r*s.sqr - h.bar/4)/(n.samp - 1)
  ht <- p.bar*(1-p.bar) + s.sqr/r
  hb/ht
}

#return an averaged Nei's FST via ratio of averages
RoA.Fst.Nei <- function(genmat, n, ploidy=2){
  if(length(n)==1){n <- rep(n, ncol(genmat)/n)}
  af.w <- sample.af.w(genmat, n, ploidy)
  af.t <- sample.af.t(genmat, ploidy)
  ht <- af.t*(1-af.t)
  hw <- ((af.w*(1-af.w)) %*% n) / sum(n)
  (sum(ht) - sum(hw))/sum(ht)
}

#return an averaged Nei's FST via average of ratios
AoR.Fst.Nei <- function(genmat, n, ploidy=2){
  lFsts <- locus.Fsts.Nei(genmat, n, ploidy)
  mean(lFsts[!is.nan(lFsts)])
}

#return an averaged Weir-Cockerham FST via ratio of averages
RoA.Fst.WC <- function(genmat, n, ploidy=2){
  if(length(n)==1){n <- rep(n, ncol(genmat)/n)}
  n.samp <- sum(n)/length(n)
  r <- length(n)
  p.bar <- sample.af.t(genmat, 2)
  s.sqr <- apply(sample.af.w(genmat, n, 2), 1, function(x) var(x))
  h.bar <- apply(genmat, 1, function(x) sum(x==1)/sum(n))
  hb <- s.sqr - (p.bar*(1-p.bar) - (r-1)/r*s.sqr - h.bar/4)/(n.samp - 1)
  ht <- p.bar*(1-p.bar) + s.sqr/r
  sum(hb)/sum(ht)
}

#return an averaged Weir-Cockerham FST via average of ratios
AoR.Fst.WC <- function(genmat, n, ploidy=2){
  lFsts <- locus.Fsts.WC(genmat, n, ploidy)
  mean(lFsts[!is.nan(lFsts)])
}

#a helper function to compute split time or migration rate 
#for a model with specified equilibrium Nei's FST
getPar <- function(Ne, d, FST){
  if (model == "split_star") {return(FST/(1-FST) * 2*Ne*d/(d-1))}
  if (model == "split_balanced") {return(FST/(1-FST) * 2*Ne*d/(d*log2(d)-d+1))}
  if (model == "split_caterpillar") {return(FST/(1-FST) * 6*Ne*d/((d-1)*(2*d-1)))}
  if (model == "migration_island") {return((1-FST)/FST * ((d-1)**2)/(4*Ne*d**2))}
  if (model == "migration_circular") {return((1-FST)/FST * ((d**2)-1)/(24*Ne*d))}
}

#construct a coalescence matrix for a star-like split model
#of d demes, with a population size of Ne per deme and a split time t
star.matrix <- function(Ne, d, t){
  res <- matrix(2*Ne, d, d)
  for (i in 1:d) {
    for (j in 1:d) {
      if(i != j){
        res[i, j] = res[i, j] + t
      }
    }
  }
  return(res)
}

#return a coalescence matrix for a balanced split model
balanced.matrix <- function(Ne, d, t){
  res <- matrix(2*Ne,d,d)
  for (i in 1:log2(d)){
    res[1:(2^(i-1)), (2^(i-1)+1):(2^i)] = res[1:(2^(i-1)), (2^(i-1)+1):(2^i)] + i*t
    res[(2^(i-1)+1):(2^i), 1:(2^(i-1))] = res[(2^(i-1)+1):(2^i), 1:(2^(i-1))] + i*t
    res[(2^(i-1)+1):(2^i), (2^(i-1)+1):(2^i)] = res[1:(2^(i-1)), 1:(2^(i-1))]
  }
  return(res)
}

#return a coalescence matrix for a caterpillar split model
caterpillar.matrix <- function(Ne, d, t){
  res <- matrix(2*Ne, d, d)
  for (i in 1:d) {
    for (j in 1:d) {
      if(i != j){
        res[i, j] = res[i, j] + (d-min(i,j))*t
      }
    }
  }
  return(res)
}

#given the number of demes, effective population size, and migration rate
#return a coalescence matrix for an island model
island.matrix <- function(Ne, d, m){
  res <- matrix(2*Ne*d, d, d)
  for (i in 1:d) {
    for (j in 1:d) {
      if(i != j){
        res[i, j] = res[i, j] + (d-1)/(2*m)
      }
    }
  }
  return(res)
}

#return a coalescence matrix for a circular stepping-stone model
circular.matrix <- function(Ne, d, m){
  res <- matrix(2*Ne*d, d, d)
  for (i in 1:d) {
    for (j in 1:d) {
      b <- abs(i - j)
      res[i,j] = res[i,j] + (d-b)*b/(2*m)
    }
  }
  return(res)
}

#a helper function to get the coalescence matrix of a specific demographic model
getMatrix <- function(Ne, d, par){
  if (model == "split_star") {return(star.matrix(Ne, d, par))}
  if (model == "split_balanced") {return(balanced.matrix(Ne, d, par))}
  if (model == "split_caterpillar") {return(caterpillar.matrix(Ne, d, par))}
  if (model == "migration_island") {return(island.matrix(Ne, d, par))}
  if (model == "migration_circular") {return(circular.matrix(Ne, d, par))}
}

#simulate a multivariate normal distribution as a null distribution for Qst
#according to Evan Koch's approach (Evan Koch, 2019)
#takes in a coalescence matrix
qst.sim <- function(coal.mat, nn=1){
  cov.mat <- coal.mat
  for(ii in 1:ncol(coal.mat))
    for(jj in 1:ncol(coal.mat))
      cov.mat[ii, jj] <- mean(coal.mat[ii,]) + mean(coal.mat[,jj]) - mean(coal.mat) - coal.mat[ii,jj]
  sq.devs <- mvrnorm(nn, rep(0, ncol(coal.mat)), Sigma=cov.mat)^2
  var.between <- apply(sq.devs, 1, mean)
  return(var.between / (var.between + mean(diag(coal.mat))))
}



######################################
# Initialization
######################################
setwd("/Users/BidYD/Desktop/F/fst")
#install.packages("vcd")
#install.packages("gridExtra")
library(vcd)
library(gridExtra)
require(MASS)
nulls <- c('single-locus Nei', 'single-locus WC', 'LK Nei RoA', 'LK Nei AoR', 'LK WC RoA', 'LK WC AoR', 'Koch MVN', 'common variants Nei', 'common variants WC')
pal <- c('#CC6677', '#117733', '#88CCEE','#882255','#999933', '#44AA99', '#332288','#DDCC77','#AA4499')




######################################
# Figure 2 & Supp: QST behavior
######################################
data <- matrix(NA, 0, 8)
for (model in c("split_star", "split_balanced", "split_caterpillar", "migration_island", "migration_circular")){
  print(model)
  for (demes in c(2,4,8,16)){
    if (any(file.exists(paste0("output/gaussian/", demes, "D_", model, "_0.1_qst_mtx.rds")))){
      print(paste0(demes, " demes"))
      Qstmat <- readRDS(file.path(paste0("output/gaussian/", demes, "D_", model, "_0.1_qst_mtx.rds")))
      data <- rbind(data, colMeans(Qstmat))
    }
  }
}
data.rb <- data[,c(1,3,5,7)]
data.pbs <- data[,c(2,4,6,8)]
rownames(data.rb) <- c("2D starlike", "4D starlike", "8D starlike", "16D starlike",
                       "4D balanced", "8D balanced", "16D balanced",
                       "4D caterpillar", "8D caterpillar", "16D caterpillar",
                       "2D island", "4D island", "8D island", "16D island",
                       "4D circular", "8D circular", "16D circular")
rownames(data.pbs) <- rownames(data.rb)


png("figures/raw/fig2/Supp_QstRB_behavior_split.png", width = 2100, height = 2100, res = 300)

par(mar=c(5,6,2,2))
plot(x = log10(c(1, 10, 100, 1000)), y = c(0,0,0,0), ylim = c(0, 0.15), pch = "", 
     bty = "n", ylab = "", xlab = "", xaxt = "n", yaxt = "n")
axis(1, at = c(0, 1, 2, 3), labels = c("1", "10", "100", "1000"), cex.axis = 2)
axis(2, las = 1, cex.axis = 2)
lines(log10(c(1, 10, 100, 1000)), data.rb[1,], type = "o", pch = 16, col = "#a8ddb5", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.rb[2,], type = "o", pch = 16, col = "#7bccc4", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.rb[3,], type = "o", pch = 16, col = "#43a2ca", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.rb[4,], type = "o", pch = 16, col = "#0868ac", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.rb[5,], type = "o", pch = 16, col = "#78c679", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.rb[6,], type = "o", pch = 16, col = "#31a354", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.rb[7,], type = "o", pch = 16, col = "#006837", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.rb[8,], type = "o", pch = 16, col = "#8c96c6", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.rb[9,], type = "o", pch = 16, col = "#8856a7", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.rb[10,], type = "o", pch = 16, col = "#810f7c", lwd = 2)
lines(c(0,100), c(.1, .1), lty = 2, col = "#bdbdbd")
title(ylab=expression(Q[ST]^RB), line=3, cex.lab=2.5)
title(xlab = "Number of causal loci", line=3, cex.lab=2.5)
legend("bottomright", bty = "n", lty = c(2, rep(1, 10)), pch = c(NA, rep(16, 10)), lwd = 2,
       col = c("#bdbdbd","#a8ddb5",
               "#7bccc4", "#78c679", "#8c96c6",
               "#43a2ca", "#31a354", "#8856a7",
               "#0868ac", "#006837", "#810f7c"),
       legend = c(expression(F[ST]^Nei ~ (t-t[W])/t),"2D",
                  "4D star-like", "4D balanced", "4D caterpillar",
                  "8D star-like", "8D balanced", "8D caterpillar",
                  "16D star-like", "16D balanced", "16D caterpillar"), cex = 1.2)
# legend("bottomright", bty = "n", lty = c(2, rep(1, 4)), pch = c(NA, rep(16, 4)), lwd = 2,
#        col = c("#bdbdbd","#a8ddb5", "#7bccc4", "#43a2ca", "#0868ac"),
#        legend = c(expression(F[ST]^Nei ~ (t-t[W])/t),"2D", "4D star-like", "8D star-like", "16D star-like"), cex = 1.5)

dev.off()


png("figures/raw/fig2/Supp_QstPBS_behavior_split.png", width = 2100, height = 2100, res = 300)

par(mar=c(5,6,2,2))
plot(x = log10(c(1, 10, 100, 1000)), y = c(0,0,0,0), ylim = c(0, 0.15), pch = "", 
     bty = "n", ylab = "", xlab = "", xaxt = "n", yaxt = "n")
axis(1, at = c(0, 1, 2, 3), labels = c("1", "10", "100", "1000"), cex.axis = 2)
axis(2, las = 1, cex.axis = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs[1,], type = "o", pch = 16, col = "#a8ddb5", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs[2,], type = "o", pch = 16, col = "#7bccc4", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs[3,], type = "o", pch = 16, col = "#43a2ca", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs[4,], type = "o", pch = 16, col = "#0868ac", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs[5,], type = "o", pch = 16, col = "#78c679", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs[6,], type = "o", pch = 16, col = "#31a354", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs[7,], type = "o", pch = 16, col = "#006837", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs[8,], type = "o", pch = 16, col = "#8c96c6", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs[9,], type = "o", pch = 16, col = "#8856a7", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs[10,], type = "o", pch = 16, col = "#810f7c", lwd = 2)
lines(c(0,100), c(.1, .1), lty = 2, col = "#bdbdbd")
title(ylab=expression(Q[ST]^PBS), line=3, cex.lab=2.5)
title(xlab = "Number of causal loci", line=3, cex.lab=2.5)
legend("bottomright", bty = "n", lty = c(2, rep(1, 10)), pch = c(NA, rep(16, 10)), lwd = 2,
       col = c("#bdbdbd","#a8ddb5",
               "#7bccc4", "#78c679", "#8c96c6",
               "#43a2ca", "#31a354", "#8856a7",
               "#0868ac", "#006837", "#810f7c"),
       legend = c(expression(F[ST]^Nei ~ (t-t[W])/t),"2D",
                  "4D star-like", "4D balanced", "4D caterpillar",
                  "8D star-like", "8D balanced", "8D caterpillar",
                  "16D star-like", "16D balanced", "16D caterpillar"), cex = 1.2)
# legend("bottomright", bty = "n", lty = c(2, rep(1, 4)), pch = c(NA, rep(16, 4)), lwd = 2,
#        col = c("#bdbdbd","#a8ddb5", "#7bccc4", "#43a2ca", "#0868ac"),
#        legend = c(expression(F[ST]^Nei ~ (t-t[W])/t),"2D", "4D star-like", "8D star-like", "16D star-like"), cex = 1.5)

dev.off()


png("figures/raw/fig2/Supp_QstRB_behavior_migration.png", width = 2100, height = 2100, res = 300)

par(mar=c(5,6,2,2))
plot(x = log10(c(1, 10, 100, 1000)), y = c(0,0,0,0), ylim = c(0, 0.15), pch = "", 
     bty = "n", ylab = "", xlab = "", xaxt = "n", yaxt = "n")
axis(1, at = c(0, 1, 2, 3), labels = c("1", "10", "100", "1000"), cex.axis = 2)
axis(2, las = 1, cex.axis = 2)
lines(log10(c(1, 10, 100, 1000)), data.rb[11,], type = "o", pch = 16, col = "#fec44f", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.rb[12,], type = "o", pch = 16, col = "#fe9929", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.rb[13,], type = "o", pch = 16, col = "#d95f0e", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.rb[14,], type = "o", pch = 16, col = "#993404", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.rb[15,], type = "o", pch = 16, col = "#f768a1", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.rb[16,], type = "o", pch = 16, col = "#c51b8a", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.rb[17,], type = "o", pch = 16, col = "#7a0177", lwd = 2)
lines(c(0,100), c(.1, .1), lty = 2, col = "#bdbdbd")
title(ylab=expression(Q[ST]^RB), line=3, cex.lab=2.5)
title(xlab = "Number of causal loci", line=3, cex.lab=2.5)
legend("bottomright", bty = "n", lty = c(2, rep(1, 7)), pch = c(NA, rep(16, 7)), lwd = 2,
       col = c("#bdbdbd", "#fec44f",
               "#fe9929", "#f768a1",
               "#d95f0e", "#c51b8a",
               "#993404", "#7a0177"),
       legend = c(expression(F[ST]^Nei ~ (t-t[W])/t),"2D",
                  "4D island", "4D circular",
                  "8D island", "8D circular",
                  "16D island", "16D circular"), cex = 1.2)
# legend("bottomright", bty = "n", lty = c(2, rep(1, 4)), pch = c(NA, rep(16, 4)), lwd = 2,
#        col = c("#bdbdbd", "#fec44f", "#fe9929", "#d95f0e", "#993404"),
#        legend = c(expression(F[ST]^Nei ~ (t-t[W])/t),"2D", "4D island", "8D island", "16D island"), cex = 1.5)

dev.off()



png("figures/raw/fig2/Supp_QstPBS_behavior_migration.png", width = 2100, height = 2100, res = 300)

par(mar=c(5,6,2,2))
plot(x = log10(c(1, 10, 100, 1000)), y = c(0,0,0,0), ylim = c(0, 0.15), pch = "", 
     bty = "n", ylab = "", xlab = "", xaxt = "n", yaxt = "n")
axis(1, at = c(0, 1, 2, 3), labels = c("1", "10", "100", "1000"), cex.axis = 2)
axis(2, las = 1, cex.axis = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs[11,], type = "o", pch = 16, col = "#fec44f", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs[12,], type = "o", pch = 16, col = "#fe9929", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs[13,], type = "o", pch = 16, col = "#d95f0e", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs[14,], type = "o", pch = 16, col = "#993404", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs[15,], type = "o", pch = 16, col = "#f768a1", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs[16,], type = "o", pch = 16, col = "#c51b8a", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs[17,], type = "o", pch = 16, col = "#7a0177", lwd = 2)
lines(c(0,100), c(.1, .1), lty = 2, col = "#bdbdbd")
title(ylab=expression(Q[ST]^PBS), line=3, cex.lab=2.5)
title(xlab = "Number of causal loci", line=3, cex.lab=2.5)
legend("bottomright", bty = "n", lty = c(2, rep(1, 7)), pch = c(NA, rep(16, 7)), lwd = 2,
       col = c("#bdbdbd", "#fec44f",
               "#fe9929", "#f768a1",
               "#d95f0e", "#c51b8a",
               "#993404", "#7a0177"),
       legend = c(expression(F[ST]^Nei ~ (t-t[W])/t),"2D",
                  "4D island", "4D circular",
                  "8D island", "8D circular",
                  "16D island", "16D circular"), cex = 1.2)
# legend("bottomright", bty = "n", lty = c(2, rep(1, 4)), pch = c(NA, rep(16, 4)), lwd = 2,
#        col = c("#bdbdbd", "#fec44f", "#fe9929", "#d95f0e", "#993404"),
#        legend = c(expression(F[ST]^Nei ~ (t-t[W])/t),"2D", "4D island", "8D island", "16D island"), cex = 1.5)

dev.off()


FST.WC.2 <- 2/11
FST.WC.4 <- 4/31
FST.WC.8 <- 8/71
FST.WC.16 <- 16/151
data.pbs.r <- cbind(data.pbs, c(FST.WC.2, FST.WC.4, FST.WC.8, FST.WC.16,
                                FST.WC.4, FST.WC.8, FST.WC.16,
                                FST.WC.4, FST.WC.8, FST.WC.16,
                                FST.WC.2, FST.WC.4, FST.WC.8, FST.WC.16,
                                FST.WC.4, FST.WC.8, FST.WC.16))
data.pbs.r <- data.pbs.r/data.pbs.r[,5]
data.pbs.r <- data.pbs.r[,1:4]


png("figures/raw/fig2/Supp_QstPBSratio_behavior_split.png", width = 2100, height = 2100, res = 300)

par(mar=c(5,7,2,2))
plot(x = log10(c(1, 10, 100, 1000)), y = c(0,0,0,0), ylim = c(0, 1.2), pch = "", 
     bty = "n", ylab = "", xlab = "", xaxt = "n", yaxt = "n")
axis(1, at = c(0, 1, 2, 3), labels = c("1", "10", "100", "1000"), cex.axis = 2)
axis(2, las = 1, cex.axis = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs.r[1,], type = "o", pch = 16, col = "#a8ddb5", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs.r[2,], type = "o", pch = 16, col = "#7bccc4", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs.r[3,], type = "o", pch = 16, col = "#43a2ca", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs.r[4,], type = "o", pch = 16, col = "#0868ac", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs.r[5,], type = "o", pch = 16, col = "#78c679", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs.r[6,], type = "o", pch = 16, col = "#31a354", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs.r[7,], type = "o", pch = 16, col = "#006837", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs.r[8,], type = "o", pch = 16, col = "#8c96c6", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs.r[9,], type = "o", pch = 16, col = "#8856a7", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs.r[10,], type = "o", pch = 16, col = "#810f7c", lwd = 2)
lines(c(0,100), c(1, 1), lty = 2, col = "#bdbdbd")
title(ylab=expression(Q[ST]^PBS/F[ST]^WC), line=4, cex.lab=2.5)
title(xlab = "Number of causal loci", line=3, cex.lab=2.5)
legend("bottomright", bty = "n", lty = c(2, rep(1, 10)), pch = c(NA, rep(16, 10)), lwd = 2,
       col = c("#bdbdbd","#a8ddb5",
               "#7bccc4", "#78c679", "#8c96c6",
               "#43a2ca", "#31a354", "#8856a7",
               "#0868ac", "#006837", "#810f7c"),
       legend = c(expression("expected"~Q[ST]^PBS/F[ST]^WC),"2D",
                  "4D star-like", "4D balanced", "4D caterpillar",
                  "8D star-like", "8D balanced", "8D caterpillar",
                  "16D star-like", "16D balanced", "16D caterpillar"), cex = 1.2)

dev.off()


png("figures/raw/fig2/Supp_QstPBSratio_behavior_migration.png", width = 2100, height = 2100, res = 300)

par(mar=c(5,7,2,2))
plot(x = log10(c(1, 10, 100, 1000)), y = c(0,0,0,0), ylim = c(0, 1.2), pch = "", 
     bty = "n", ylab = "", xlab = "", xaxt = "n", yaxt = "n")
axis(1, at = c(0, 1, 2, 3), labels = c("1", "10", "100", "1000"), cex.axis = 2)
axis(2, las = 1, cex.axis = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs.r[11,], type = "o", pch = 16, col = "#fec44f", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs.r[12,], type = "o", pch = 16, col = "#fe9929", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs.r[13,], type = "o", pch = 16, col = "#d95f0e", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs.r[14,], type = "o", pch = 16, col = "#993404", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs.r[15,], type = "o", pch = 16, col = "#f768a1", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs.r[16,], type = "o", pch = 16, col = "#c51b8a", lwd = 2)
lines(log10(c(1, 10, 100, 1000)), data.pbs.r[17,], type = "o", pch = 16, col = "#7a0177", lwd = 2)
lines(c(0,100), c(1, 1), lty = 2, col = "#bdbdbd")
title(ylab=expression(Q[ST]^PBS/F[ST]^WC), line=4, cex.lab=2.5)
title(xlab = "Number of causal loci", line=3, cex.lab=2.5)
legend("bottomright", bty = "n", lty = c(2, rep(1, 7)), pch = c(NA, rep(16, 7)), lwd = 2,
       col = c("#bdbdbd", "#fec44f",
               "#fe9929", "#f768a1",
               "#d95f0e", "#c51b8a",
               "#993404", "#7a0177"),
       legend = c(expression("expected"~Q[ST]^PBS/F[ST]^WC),"2D",
                  "4D island", "4D circular",
                  "8D island", "8D circular",
                  "16D island", "16D circular"), cex = 1.2)

dev.off()




######################################
# Figure 3: Single-locus FST nulls
######################################
gens <- data.matrix(read.table(paste0("input/fst0.1/geno_8D_migration_island.txt"), sep = ""))
Qstmat <- readRDS(file.path(paste0("output/gaussian/8D_migration_island_0.1_qst_mtx.rds")))
n <- rep(2000, 4)
Fsts.Nei <- locus.Fsts.Nei(gens, n)
Fsts.WC <- locus.Fsts.WC(gens, n)


pdf(paste0("figures/Figure_3.pdf"), width = 7, height = 7)

par(oma = c(4,4,1,1), mfrow=c(2,2), mar=c(2,2,1,1))
for (i in 1:4){
  hist(Qstmat[,i*2-1], freq = FALSE, breaks = 100, main = "", bty = "n",
       ylim=c(0,35),
       xlab = "", ylab = "", yaxt = "n", cex.lab=2, cex.axis=1.5)
  axis(2, las = 1, cex.axis = 1.5)
  lines(density(Fsts.Nei), col = '#CC6677', lwd = 2)
  lines(density(Fsts.WC), col = '#117733', lwd = 2)
  l = 10^(i-1)
  if (l==1) lable="locus" else lable="loci"
  legend("topright", bty = "n", lwd = 0,
         legend = bquote(.(l) ~ .(lable)), cex = 1.2)
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
title(ylab="Density", line=-2, cex.lab=2)
title(xlab = expression(Q[ST]^RB), line=-3, cex.lab=2)
legend("bottom", bty = "n", col = c('#CC6677', '#117733'), lwd = 2, horiz = TRUE, xpd = TRUE,
       legend = c("single-locus Nei", "single-locus WC"),
       cex = 1.5)

dev.off()



######################################
# Figure 4: LK nulls
######################################
gens <- data.matrix(read.table(paste0("input/fst0.1/geno_8D_split_star.txt"), sep = ""))
Qstmat <- readRDS(file.path(paste0("output/gaussian/8D_split_star_0.1_qst_mtx.rds")))
n <- rep(2000, 8)
Nei.RoA <- RoA.Fst.Nei(gens, n)
Nei.AoR <- AoR.Fst.Nei(gens, n)
WC.RoA <- RoA.Fst.WC(gens, n)
WC.AoR <- AoR.Fst.WC(gens, n)
Fstat <- seq(0,1, by = .0001)
lknra.dens <- dchisq((8 - 1) * Fstat / Nei.RoA, df = 8 - 1) * (8 - 1) / Nei.RoA
lknar.dens <- dchisq((8 - 1) * Fstat / Nei.AoR, df = 8 - 1) * (8 - 1) / Nei.AoR
lkwra.dens <- dchisq((8 - 1) * Fstat / WC.RoA, df = 8 - 1) * (8 - 1) / WC.RoA
lkwar.dens <- dchisq((8 - 1) * Fstat / WC.AoR, df = 8 - 1) * (8 - 1) / WC.AoR


pdf(paste0("figures/Figure_4.pdf"), width = 7, height = 7)

par(oma = c(4,4,1,1), mfrow=c(2,2), mar=c(2,2,1,1))
for (i in 1:4){
  hist(Qstmat[,i*2-1], freq = FALSE, breaks = 100, main = "", bty = "n",
       ylim=c(0,35),
       xlab = "", ylab = "", yaxt = "n", cex.lab=2, cex.axis=1.5)
  axis(2, las = 1, cex.axis = 1.5)
  lines(Fstat, lknra.dens, col = '#88CCEE', lwd = 2)
  lines(Fstat, lknar.dens, col = '#882255', lwd = 2)
  lines(Fstat, lkwra.dens, col = '#999933', lwd = 2)
  lines(Fstat, lkwar.dens, col = '#44AA99', lwd = 2)
  l = 10^(i-1)
  if (l==1) lable="locus" else lable="loci"
  legend("topright", bty = "n", lwd = 0,
         legend = bquote(.(l) ~ .(lable)), cex = 1.2)
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
title(ylab="Density", line=-2, cex.lab=2)
title(xlab = expression(Q[ST]^RB), line=-3, cex.lab=2)
legend("bottom", bty = "n", col = c('#88CCEE','#882255','#999933','#44AA99'), lwd = 2, horiz = TRUE, xpd = TRUE,
       legend = c('LK Nei RoA', 'LK Nei AoR', 'LK WC RoA', 'LK WC AoR'),
       cex = 1.2)

dev.off()




######################################
# Figure 5: Comparison of nulls in a 
# spatial structured demography
######################################
Fstat <- seq(0,1, by = .0001)
for (model in c("migration_circular")){
  print(model)
  for (demes in c(4,8,16)){
    if (any(file.exists(paste0("output/gaussian/", demes, "D_", model, "_0.1_qst_mtx.rds")))){
      print(paste0(demes, " demes"))
      
      gens <- data.matrix(read.table(paste0("input/fst0.1/geno_", demes, "D_", model, ".txt"), sep = ""))
      Qstmat <- readRDS(file.path(paste0("output/gaussian/", demes, "D_", model, "_0.1_qst_mtx.rds")))
      n <- rep(2000, demes)
      
      gens.freq <- cbind(gens, "freq" = sample.af.t(gens))
      gens.cv <- gens.freq[which(gens.freq[,ncol(gens.freq)]>0.05),-1]
      Fsts.Nei.cv <- locus.Fsts.Nei(gens.cv, n)
      
      Nei.RoA <- RoA.Fst.Nei(gens, n)
      lknra.dens <- dchisq((demes - 1) * Fstat / Nei.RoA, df = demes - 1) * (demes - 1) / Nei.RoA
      
      coal.mtx <- getMatrix(1000, demes, getPar(1000, demes, 0.1))
      Koch.mvn <- qst.sim(coal.mtx, nn=1e6)
      
      locusNeiCv.cut <- quantile(Fsts.Nei.cv, 0.95, na.rm = TRUE)
      lknra.cut <- qchisq(0.95, demes - 1) * Nei.RoA / (demes - 1)
      Mvn.cut <- quantile(Koch.mvn, 0.95, na.rm = TRUE)
      
      
      png(paste0("figures/raw/fig5/", demes, "D_", model, "_QST_RB.png"),
          width = 2100, height = 2100, res = 300)
      
      par(mar=c(5,6,2,2))
      hist(Qstmat[,7], freq = FALSE, breaks = 100, main = "", bty = "n",
           ylim=c(0,15),
           xlab = "", ylab = "", yaxt = "n", cex.lab=2, cex.axis=2)
      axis(2, las = 1, cex.axis = 2)
      lines(density(Fsts.Nei.cv), col = '#DDCC77', lwd = 2)
      lines(Fstat, lknra.dens, col = '#88CCEE', lwd = 2)
      lines(density(Koch.mvn), col = '#332288', lwd = 2)
      title(ylab="Density", line=4, cex.lab=2.5)
      title(xlab = expression(Q[ST]^RB), line=4, cex.lab=2.5)
      legend("topright", bty = "n", lwd = 0,
             legend = bquote(.(demes) ~ "demes"), cex = 1.6)
      legend("topleft", bty = "n", col = c('#DDCC77','#88CCEE','#332288'), lwd = 2,
             legend = c('common variants Nei', 'LK Nei RoA', 'Koch MVN'),
             cex = 1.6)
      
      dev.off()
      
      
      png(paste0("figures/raw/fig5/", demes, "D_", model, "_error_rates.png"),
          width = 2100, height = 2100, res = 300)
      
      par(mar=c(5,7,2,2))
      plot(x = log10(c(1, 10, 100, 1000)), y = c(0,0,0,0), ylim = c(0, 0.15), pch = "",
           bty = "n", ylab = "", xlab = "", xaxt = "n", yaxt = "n")
      axis(1, at = c(0, 1, 2, 3), labels = c("1", "10", "100", "1000"), cex.axis = 2)
      axis(2, las = 1, cex.axis = 2)
      points(log10(c(1, 10, 100, 1000)) - 0.05, colMeans(Qstmat > locusNeiCv.cut)[c(1,3,5,7)], col = '#DDCC77', pch = 19)
      points(log10(c(1, 10, 100, 1000)), colMeans(Qstmat > lknra.cut)[c(1,3,5,7)], col = '#88CCEE', pch = 19)
      points(log10(c(1, 10, 100, 1000)) + 0.05, colMeans(Qstmat > Mvn.cut)[c(1,3,5,7)], col = '#332288', pch = 19)
      lines(c(0,100), c(.05, .05), lty = 2, col = "grey")
      title(ylab="Type I error rate", line=5, cex.lab=2.5)
      title(xlab = "Number of causal loci", line=3, cex.lab=2.5)
      legend("topleft", bty = "n", col = c('#DDCC77','#88CCEE','#332288'), pch = c(19, 19, 19),
             legend = c('common variants Nei', 'LK Nei RoA', 'Koch MVN'), cex = 1.6)
      legend("topright", bty = "n", lwd = 0,
             legend = bquote(.(demes) ~ "demes"), cex = 1.6)

      
      se <- sqrt(colMeans(Qstmat > locusNeiCv.cut)*(1 - colMeans(Qstmat > locusNeiCv.cut))/nrow(Qstmat))
      lb <-colMeans(Qstmat > locusNeiCv.cut) - 2*se
      ub <-colMeans(Qstmat > locusNeiCv.cut) + 2*se
      for(i in c(1,3,5,7)){
        lines(rep(floor((i - 1)/2)  - .05, 2), c(lb[i], ub[i]), col = '#DDCC77', lwd = 2)
      }
      
      se <- sqrt(colMeans(Qstmat > lknra.cut)*(1 - colMeans(Qstmat > lknra.cut))/nrow(Qstmat))
      lb <-colMeans(Qstmat > lknra.cut) - 2*se
      ub <-colMeans(Qstmat > lknra.cut) + 2*se
      for(i in c(1,3,5,7)){
        lines(rep(floor((i - 1)/2), 2), c(lb[i], ub[i]), col = '#88CCEE', lwd = 2)
      }
      
      se <- sqrt(colMeans(Qstmat > Mvn.cut)*(1 - colMeans(Qstmat > Mvn.cut))/nrow(Qstmat))
      lb <-colMeans(Qstmat > Mvn.cut) - 2*se
      ub <-colMeans(Qstmat > Mvn.cut) + 2*se
      for(i in c(1,3,5,7)){
        lines(rep(floor((i - 1)/2) + .05, 2), c(lb[i], ub[i]), col = '#332288', lwd = 2)
      }
      
      dev.off()
    }
  }
}




######################################
# Figure 6: Error rates
######################################
### Panel A: LK Nei RoA vs AoR - QST RB
getER.A <- function(Qstvec, lknra.cut, lknar.cut){
  er1 <- colMeans(Qstvec > lknra.cut)
  er2 <- colMeans(Qstvec > lknar.cut)
  
  se1 <- sqrt(er1*(1 - er1)/nrow(Qstvec))
  lb1 <- er1 - 2*se1
  ub1 <- er1 + 2*se1
  se2 <- sqrt(er2*(1 - er2)/nrow(Qstvec))
  lb2 <- er2 - 2*se2
  ub2 <- er2 + 2*se2
  
  ub.y <- max(er1[c(1,3,5,7)], er2[c(1,3,5,7)])*5/4
  if (ub.y > 1){ub.y = 1}
  par(mar=c(5,7,2,2))
  plot(x = log10(c(1, 10, 100, 1000)), y = c(0,0,0,0), ylim = c(0, ub.y), pch = "",
       bty = "n", ylab = "", xlab = "", xaxt = "n", yaxt = "n")
  axis(1, at = c(0, 1, 2, 3), labels = c("1", "10", "100", "1000"), cex.axis = 2)
  axis(2, las = 1, cex.axis = 2)
  points(log10(c(1, 10, 100, 1000)) - 0.025, er1[c(1,3,5,7)], col = '#88CCEE', pch = 19)
  points(log10(c(1, 10, 100, 1000)) + 0.025, er2[c(1,3,5,7)], col = '#882255', pch = 19)
  lines(c(0,100), c(.05, .05), lty = 2, col = "grey")
  title(ylab="Type I error rate", line=5, cex.lab=2.5)
  title(xlab = "Number of causal loci", line=3, cex.lab=2.5)
  legend("topleft", bty = "n", col = c('#88CCEE','#882255'), pch = c(19, 19),
         legend = c('LK Nei RoA', 'LK Nei AoR'), cex = 1.6)
  legend("topright", bty = "n", lwd = 0,
         legend = bquote(.(demes) ~ "demes"), cex = 1.6)
  for(i in c(1,3,5,7)){
    lines(rep(floor((i - 1)/2) - .025, 2), c(lb1[i], ub1[i]), col = '#88CCEE', lwd = 2)
    lines(rep(floor((i - 1)/2) + .025, 2), c(lb2[i], ub2[i]), col = '#882255', lwd = 2)
  }
}


### Panel B: LK Nei RoA - QST RB vs QST PBS
getER.B <- function(Qstvec, lknra.cut){
  er <- colMeans(Qstvec > lknra.cut)
  se <- sqrt(er*(1 - er)/nrow(Qstvec))
  lb <- er - 2*se
  ub <- er + 2*se
  
  ub.y <- max(er)*5/4
  if (ub.y > 1){ub.y = 1}
  par(mar=c(5,7,2,2))
  plot(x = log10(c(1, 10, 100, 1000)), y = c(0,0,0,0), ylim = c(0, ub.y), pch = "",
       bty = "n", ylab = "", xlab = "", xaxt = "n", yaxt = "n")
  axis(1, at = c(0, 1, 2, 3), labels = c("1", "10", "100", "1000"), cex.axis = 2)
  axis(2, las = 1, cex.axis = 2)
  points(log10(c(1, 10, 100, 1000)), er[c(1,3,5,7)], col = '#88CCEE', pch = 19)
  points(log10(c(1, 10, 100, 1000)), er[c(2,4,6,8)], col = '#88CCEE', pch = 9)
  lines(c(0,100), c(.05, .05), lty = 2, col = "grey")
  title(ylab="Type I error rate", line=5, cex.lab=2.5)
  title(xlab = "Number of causal loci", line=3, cex.lab=2.5)
  legend("topleft", bty = "n", col = c('#88CCEE','black','black'), pch = c(19,19,9),
         legend = c('LK Nei RoA', expression("RB"~Q[ST]), expression("PBS"~Q[ST])), cex = 1.6)
  legend("topright", bty = "n", lwd = 0,
         legend = bquote(.(demes) ~ "demes"), cex = 1.6)
  for(i in 1:8){
    lines(rep(floor((i - 1)/2), 2), c(lb[i], ub[i]), col = '#88CCEE', lwd = 2)
  }
}


### Panel C: single-locus Nei vs common variants Nei vs Koch MVN - QST RB
getER.C <- function(Qstvec, locusNei.cut, locusNeiCv.cut, Mvn.cut){
  er1 <- colMeans(Qstvec > locusNei.cut)
  er2 <- colMeans(Qstvec > locusNeiCv.cut)
  er3 <- colMeans(Qstvec > Mvn.cut)
  
  se1 <- sqrt(er1*(1 - er1)/nrow(Qstvec))
  lb1 <- er1 - 2*se1
  ub1 <- er1 + 2*se1
  se2 <- sqrt(er2*(1 - er2)/nrow(Qstvec))
  lb2 <- er2 - 2*se2
  ub2 <- er2 + 2*se2
  se3 <- sqrt(er3*(1 - er3)/nrow(Qstvec))
  lb3 <- er3 - 2*se3
  ub3 <- er3 + 2*se3
  
  ub.y <- max(er1[c(1,3,5,7)],er2[c(1,3,5,7)],er3[c(1,3,5,7)])*5/4
  if (ub.y > 1){ub.y = 1}
  par(mar=c(5,7,2,2))
  plot(x = log10(c(1, 10, 100, 1000)), y = c(0,0,0,0), ylim = c(0, ub.y), pch = "",
       bty = "n", ylab = "", xlab = "", xaxt = "n", yaxt = "n")
  axis(1, at = c(0, 1, 2, 3), labels = c("1", "10", "100", "1000"), cex.axis = 2)
  axis(2, las = 1, cex.axis = 2)
  points(log10(c(1, 10, 100, 1000)) - 0.05, er1[c(1,3,5,7)], col = '#CC6677', pch = 19)
  points(log10(c(1, 10, 100, 1000)), er2[c(1,3,5,7)], col = '#DDCC77', pch = 19)
  points(log10(c(1, 10, 100, 1000)) + 0.05, er3[c(1,3,5,7)], col = '#332288', pch = 19)
  lines(c(0,100), c(.05, .05), lty = 2, col = "grey")
  title(ylab="Type I error rate", line=5, cex.lab=2.5)
  title(xlab = "Number of causal loci", line=3, cex.lab=2.5)
  legend("topleft", bty = "n", col = c('#CC6677','#DDCC77','#332288'), pch = c(19,19,19),
         legend = c('single-locus Nei', 'common variants Nei', 'Koch MVN'), cex = 1.6)
  legend("topright", bty = "n", lwd = 0,
         legend = bquote(.(demes) ~ "demes"), cex = 1.6)
  for(i in c(1,3,5,7)){
    lines(rep(floor((i - 1)/2) - .05, 2), c(lb1[i], ub1[i]), col = '#CC6677', lwd = 2)
    lines(rep(floor((i - 1)/2), 2), c(lb2[i], ub2[i]), col = '#DDCC77', lwd = 2)
    lines(rep(floor((i - 1)/2) + .05, 2), c(lb3[i], ub3[i]), col = '#332288', lwd = 2)
  }
}


Fstat <- seq(0,1, by = .0001)
for (model in c("split_star")){
  print(model)
  for (demes in c(2,8,16)){
    if (any(file.exists(paste0("output/gaussian/", demes, "D_", model, "_0.1_qst_mtx.rds")))){
      print(paste0(demes, " demes"))
      
      gens <- data.matrix(read.table(paste0("input/fst0.1/geno_", demes, "D_", model, ".txt"), sep = ""))
      Qstmat <- readRDS(file.path(paste0("output/gaussian/", demes, "D_", model, "_0.1_qst_mtx.rds")))
      n <- rep(2000, demes)
      
      Nei.RoA <- RoA.Fst.Nei(gens, n)
      Nei.AoR <- AoR.Fst.Nei(gens, n)
      Fsts.Nei <- locus.Fsts.Nei(gens, n)
      gens.freq <- cbind(gens, "freq" = sample.af.t(gens))
      gens.cv <- gens.freq[which(gens.freq[,ncol(gens.freq)]>0.05),-1]
      Fsts.Nei.cv <- locus.Fsts.Nei(gens.cv, n)
      coal.mtx <- getMatrix(1000, demes, getPar(1000, demes, 0.1))
      Qst.mvn <- qst.sim(coal.mtx, nn=1e6)
      
      lknra.cut <- qchisq(0.95, demes - 1) * Nei.RoA / (demes - 1)
      lknar.cut <- qchisq(0.95, demes - 1) * Nei.AoR / (demes - 1)
      locusNei.cut <- quantile(Fsts.Nei, 0.95, na.rm = TRUE)
      locusNeiCv.cut <- quantile(Fsts.Nei.cv, 0.95, na.rm = TRUE)
      Mvn.cut <- quantile(Qst.mvn, 0.95, na.rm = TRUE)
      
      png(paste0("figures/raw/fig6/", demes, "D_", model, "_error_rates_A.png"),
          width = 2100, height = 2100, res = 300)
      getER.A(Qstmat, lknra.cut, lknar.cut)
      dev.off()
      
      png(paste0("figures/raw/fig6/", demes, "D_", model, "_error_rates_B.png"),
          width = 2100, height = 2100, res = 300)
      getER.B(Qstmat, lknra.cut)
      dev.off()
      
      png(paste0("figures/raw/fig6/", demes, "D_", model, "_error_rates_C.png"),
          width = 2100, height = 2100, res = 300)
      getER.C(Qstmat, locusNei.cut, locusNeiCv.cut, Mvn.cut)
      dev.off()
    }
  }
}




######################################
# Supp Figure: FST statistics check
######################################
getFSTdt <- function(fst){
  data.fst <- matrix(NA, 0, 12)
  for (model in c("split_star", "split_balanced", "split_caterpillar", "migration_island", "migration_circular")){
    print(model)
    for (demes in c(2,4,8,16)){
      if (any(file.exists(paste0("input/fst", fst, "/geno_", demes, "D_", model, ".txt")))){
        print(paste0(demes, " demes"))
        gens <- data.matrix(read.table(paste0("input/fst", fst, "/geno_", demes, "D_", model, ".txt"), sep = ""))
        gens.freq <- cbind(gens, "freq" = sample.af.t(gens))
        gens.ab <- gens.freq
        gens.ab[,1:2000][gens.ab[,ncol(gens.ab)] < 0.05] <- 0
        gens.ab <- gens.ab[,-1]
        gens.cv <- gens.freq[which(gens.freq[,ncol(gens.freq)] > 0.05),-1]
        n <- rep(2000, demes)
        
        Nei.AoR <- AoR.Fst.Nei(gens, n)
        Nei.RoA <- RoA.Fst.Nei(gens, n)
        WC.AoR <- AoR.Fst.WC(gens, n)
        WC.RoA <- RoA.Fst.WC(gens, n)
        AB.Nei.AoR <- AoR.Fst.Nei(gens.ab, n)
        AB.Nei.RoA <- RoA.Fst.Nei(gens.ab, n)
        AB.WC.AoR <- AoR.Fst.WC(gens.ab, n)
        AB.WC.RoA <- RoA.Fst.WC(gens.ab, n)
        CV.Nei.AoR <- AoR.Fst.Nei(gens.cv, n)
        CV.Nei.RoA <- RoA.Fst.Nei(gens.cv, n)
        CV.WC.AoR <- AoR.Fst.WC(gens.cv, n)
        CV.WC.RoA <- RoA.Fst.WC(gens.cv, n)
        
        data.fst <- rbind(data.fst, c(Nei.AoR, Nei.RoA, WC.AoR, WC.RoA, 
                                      AB.Nei.AoR, AB.Nei.RoA, AB.WC.AoR, AB.WC.RoA,
                                      CV.Nei.AoR, CV.Nei.RoA, CV.WC.AoR, CV.WC.RoA))
      }
    }
  }
  rownames(data.fst) <- c("2D starlike", "4D starlike", "8D starlike", "16D starlike",
                          "4D balanced", "8D balanced", "16D balanced",
                          "4D caterpillar", "8D caterpillar", "16D caterpillar",
                          "2D island", "4D island", "8D island", "16D island",
                          "4D circular", "8D circular", "16D circular")
  return(t(data.fst))
}

data.fst.02 <- getFSTdt(fst = 0.02)
saveRDS(data.fst.02, file.path(paste0("figures/raw/fst_stats/fst0.02_summary_statistics.rds")))
data.fst.1 <- getFSTdt(fst = 0.1)
saveRDS(data.fst.1, file.path(paste0("figures/raw/fst_stats/fst0.1_summary_statistics.rds")))
data.fst.25 <- getFSTdt(fst = 0.25)
saveRDS(data.fst.25, file.path(paste0("figures/raw/fst_stats/fst0.25_summary_statistics.rds")))

# data.fst.02 <- readRDS(file.path(paste0("figures/raw/fst_stats/fst0.02_summary_statistics.rds")))
# data.fst.1 <- readRDS(file.path(paste0("figures/raw/fst_stats/fst0.1_summary_statistics.rds")))
# data.fst.25 <- readRDS(file.path(paste0("figures/raw/fst_stats/fst0.25_summary_statistics.rds")))


getFSTsum <- function(data.fst, fst){
  for (i in c(1,5,9)){
    data.fst.sl <- data.fst[i:(i+3),]
    if (i == 1){
      png(paste0("figures/raw/fst_stats/fst", fst, "_summary.png"),
          width = 2100, height = 2100, res = 300)
    }
    if (i == 5){
      png(paste0("figures/raw/fst_stats/fst", fst, "_summary_ab.png"),
          width = 2100, height = 2100, res = 300)
    }
    if (i == 9){
      png(paste0("figures/raw/fst_stats/fst", fst, "_summary_cv.png"),
          width = 2100, height = 2100, res = 300)
    }
    
    par(mar=c(9,3,1,1))
    barplot(data.fst.sl, col=c('#7b3294', "#c2a5cf", "#008837", "#a6dba0"), beside=T, ylab = "", 
            ylim=c(0,2*fst), las=3, cex.axis = 1.2, cex.names = 1.2)
    lines(c(0,100), c(fst, fst), lty = 2, col = "#bdbdbd")
    g = grid_legend("topright", pch=15, frame=F, col=c('#7b3294', "#c2a5cf", "#008837", "#a6dba0"), 
                    labels=c("Nei AoR", "Nei RoA", "WC AoR", "WC RoA"), draw = F)
    grid.draw(grobTree(g, vp=viewport(x=0.85, y=0.89, angle=90)))
    
    dev.off()
  }
}

getFSTsum(data.fst.02, 0.02)
getFSTsum(data.fst.1, 0.1)
getFSTsum(data.fst.25, 0.25)




######################################
# Supp Figure: Effect size checks
######################################
for (model in c("split_star","migration_island")){
  qst.gs <- readRDS(paste0("output/gaussian/8D_", model, "_0.1_qst_mtx.rds"))
  qst.un <- readRDS(paste0("output/uniform/8D_", model, "_0.1_qst_mtx.rds"))
  qst.lp <- readRDS(paste0("output/laplace/8D_", model, "_0.1_qst_mtx.rds"))
  qst.fq <- readRDS(paste0("output/alpha/8D_", model, "_0.1_qst_mtx.rds"))
  
  qst.rb.1000 <- as.data.frame(cbind(qst.gs[,7], qst.un[,7], qst.lp[,7], qst.fq[,7]))
  qst.pbs.1000 <- as.data.frame(cbind(qst.gs[,8], qst.un[,8], qst.lp[,8], qst.fq[,8]))
  colnames(qst.rb.1000) <- c("gaussian", "uniform", "laplace", "freq")
  colnames(qst.pbs.1000) <- colnames(qst.rb.1000)
  
  
  png(paste0("figures/raw/qst_es/8D_", model, "_effect_sizes_rb.png"), 
      width = 2100, height = 2100, res = 300)
  
  plot(density(qst.rb.1000$gaussian),  main = "", bty = "n", xlab = "", ylab = "", 
       yaxt = "n", ylim=c(0,15), cex.lab=2, cex.axis=1.5, col = "#ca0020", lwd = 2)
  axis(2, las = 1, cex.axis = 1.5)
  lines(density(qst.rb.1000$uniform), col = "#f4a582", lwd = 2)
  lines(density(qst.rb.1000$laplace), col = "#92c5de", lwd = 2)
  lines(density(qst.rb.1000$freq), col = "#0571b0", lwd = 2)
  legend("topright", bty = "n", lwd = 2, cex = 1.5,
         col = c("#ca0020","#f4a582","#92c5de","#0571b0"), 
         legend = c("Gaussian", "uniform", "laplace", expression(alpha ~ "model (" ~ alpha ~ "= -1)")))
  title(ylab="Density", line=2.9, cex.lab=1.6)
  title(xlab = expression(Q[ST]^RB ~"with 1000 loci"), line=2.8, cex.lab=1.6)
  
  dev.off()
  
  
  png(paste0("figures/raw/qst_es/8D_", model, "_effect_sizes_pbs.png"), 
      width = 2100, height = 2100, res = 300)
  
  plot(density(qst.pbs.1000$gaussian),  main = "", bty = "n", xlab = "", ylab = "", 
       yaxt = "n", ylim=c(0,15), cex.lab=2, cex.axis=1.5, col = "#ca0020", lwd = 2)
  axis(2, las = 1, cex.axis = 1.5)
  lines(density(qst.pbs.1000$uniform), col = "#f4a582", lwd = 2)
  lines(density(qst.pbs.1000$laplace), col = "#92c5de", lwd = 2)
  lines(density(qst.pbs.1000$freq), col = "#0571b0", lwd = 2)
  legend("topright", bty = "n", lwd = 2, cex = 1.5,
         col = c("#ca0020","#f4a582","#92c5de","#0571b0"), 
         legend = c("Gaussian", "uniform", "laplace", expression(alpha ~ "model (" ~ alpha ~ "= -1)")))
  title(ylab="Density", line=2.9, cex.lab=1.6)
  title(xlab = expression(Q[ST]^PBS ~"with 1000 loci"), line=2.8, cex.lab=1.6)
  
  dev.off()
}




######################################
# Supp Figure: QST nulls - group 1
######################################
Fstat <- seq(0,1, by = .0001)

getLocusNull <- function(i, model){
  hist(Qstmat[,i], freq = FALSE, breaks = 50, main = "", bty = "n",
       ylim=c(0,15),
       xlab = "", ylab = "", yaxt = "n", cex.lab=2, cex.axis=1.5)
  axis(2, las = 1, cex.axis = 1.5)
  lines(density(Fsts.Nei), col = '#CC6677', lwd = 2)
  lines(density(Fsts.WC), col = '#117733', lwd = 2)
  lines(density(Fsts.Nei.cv), col = '#DDCC77', lwd = 2)
  lines(density(Fsts.WC.cv), col = '#AA4499', lwd = 2)
  title(ylab="Density", line=2.9, cex.lab=1.6)
  if (i == 1){title(xlab = expression(Q[ST]^RB ~"with 1 locus"), line=3.5, cex.lab=1.6)}
  else if (i == 2){title(xlab = expression(Q[ST]^PBS ~"with 1 locus"), line=3.5, cex.lab=1.6)}
  else if (i == 7){title(xlab = expression(Q[ST]^RB ~"with 1000 loci"), line=3.5, cex.lab=1.6)}
  else if (i == 8){title(xlab = expression(Q[ST]^PBS ~"with 1000 loci"), line=3.5, cex.lab=1.6)}
  
  if (model == "split_star"){legend("topright", bty = "n", lty = 0, legend = "star-like", cex = 1.5)}
  else if (model == "migration_island"){legend("topright", bty = "n", lty = 0, legend = "island", cex = 1.5)}
  else if (model == "migration_circular"){legend("topright", bty = "n", lty = 0, legend = "circular", cex = 1.5)}
}


demes = 2
pdf(paste0("figures/Figure_nulls_1_", demes, "D.pdf"),
    width = 10, height = 5)

par(oma = c(2,1,1,1), mfrow=c(2,4), mar=c(6,6,1,1))
for (model in c("split_star", "migration_island")){
  print(model)
  if (any(file.exists(paste0("output/gaussian/", demes, "D_", model, "_0.1_qst_mtx.rds")))){
    gens <- data.matrix(read.table(paste0("input/fst0.1/geno_", demes, "D_", model, ".txt"), sep = ""))
    Qstmat <- readRDS(file.path(paste0("output/gaussian/", demes, "D_", model, "_0.1_qst_mtx.rds")))
    n <- rep(2000, demes)
    Fsts.Nei <- locus.Fsts.Nei(gens, n)
    Fsts.WC <- locus.Fsts.WC(gens, n)
    gens.freq <- cbind(gens, "freq" = sample.af.t(gens))
    gens.cv <- gens.freq[which(gens.freq[,ncol(gens.freq)]>0.05),-1]
    Fsts.Nei.cv <- locus.Fsts.Nei(gens.cv, n)
    Fsts.WC.cv <- locus.Fsts.WC(gens.cv, n)
    
    getLocusNull(1, model)
    getLocusNull(7, model)
    getLocusNull(2, model)
    getLocusNull(8, model)
  }
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", bty = "n", col = c('#CC6677', '#117733','#DDCC77','#AA4499'), lwd = 2, horiz = TRUE, xpd = TRUE,
       legend = c("single-locus Nei", "single-locus WC", 'common variants Nei', 'common variants WC'),
       cex = 1.5)

dev.off()


for (demes in c(4,8,16)){
  print(paste0(demes, " demes"))
  pdf(paste0("figures/Figure_nulls_1_", demes, "D.pdf"),
      width = 10, height = 7.5)
  
  par(oma = c(2,1,1,1), mfrow=c(3,4), mar=c(6,6,1,1))
  for (model in c("split_star", "migration_island", "migration_circular")){
    print(model)
    if (any(file.exists(paste0("output/gaussian/", demes, "D_", model, "_0.1_qst_mtx.rds")))){
      gens <- data.matrix(read.table(paste0("input/fst0.1/geno_", demes, "D_", model, ".txt"), sep = ""))
      Qstmat <- readRDS(file.path(paste0("output/gaussian/", demes, "D_", model, "_0.1_qst_mtx.rds")))
      n <- rep(2000, demes)
      Fsts.Nei <- locus.Fsts.Nei(gens, n)
      Fsts.WC <- locus.Fsts.WC(gens, n)
      gens.freq <- cbind(gens, "freq" = sample.af.t(gens))
      gens.cv <- gens.freq[which(gens.freq[,ncol(gens.freq)]>0.05),-1]
      Fsts.Nei.cv <- locus.Fsts.Nei(gens.cv, n)
      Fsts.WC.cv <- locus.Fsts.WC(gens.cv, n)
      
      getLocusNull(1, model)
      getLocusNull(7, model)
      getLocusNull(2, model)
      getLocusNull(8, model)
    }
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend("bottom", bty = "n", col = c('#CC6677', '#117733','#DDCC77','#AA4499'), lwd = 2, horiz = TRUE, xpd = TRUE,
         legend = c("single-locus Nei", "single-locus WC", 'common variants Nei', 'common variants WC'),
         cex = 1.5)
  
  dev.off()
}




######################################
# Supp Figure: QST nulls - group 2
######################################
Fstat <- seq(0,1, by = .0001)

getLKMVN <- function(i, model){
  hist(Qstmat[,i], freq = FALSE, breaks = 50, main = "", bty = "n",
       ylim=c(0,15),
       xlab = "", ylab = "", yaxt = "n", cex.lab=2, cex.axis=1.5)
  axis(2, las = 1, cex.axis = 1.5)
  lines(Fstat, lknra.dens, col = '#88CCEE', lwd = 2)
  lines(Fstat, lknar.dens, col = '#882255', lwd = 2)
  lines(Fstat, lkwra.dens, col = '#999933', lwd = 2)
  lines(Fstat, lkwar.dens, col = '#44AA99', lwd = 2)
  lines(density(Koch.mvn), col = '#332288', lwd = 2)
  title(ylab="Density", line=2.9, cex.lab=1.6)
  if (i == 1){title(xlab = expression(Q[ST]^RB ~"with 1 locus"), line=3.5, cex.lab=1.6)}
  else if (i == 2){title(xlab = expression(Q[ST]^PBS ~"with 1 locus"), line=3.5, cex.lab=1.6)}
  else if (i == 7){title(xlab = expression(Q[ST]^RB ~"with 1000 loci"), line=3.5, cex.lab=1.6)}
  else if (i == 8){title(xlab = expression(Q[ST]^PBS ~"with 1000 loci"), line=3.5, cex.lab=1.6)}
  
  if (model == "split_star"){legend("topright", bty = "n", lty = 0, legend = "star-like", cex = 1.5)}
  else if (model == "migration_island"){legend("topright", bty = "n", lty = 0, legend = "island", cex = 1.5)}
  else if (model == "migration_circular"){legend("topright", bty = "n", lty = 0, legend = "circular", cex = 1.5)}
}


demes = 2
pdf(paste0("figures/Figure_nulls_2_", demes, "D.pdf"),
    width = 10, height = 5)

par(oma = c(2,1,1,1), mfrow=c(2,4), mar=c(6,6,1,1))
for (model in c("split_star", "migration_island")){
  print(model)
  if (any(file.exists(paste0("output/gaussian/", demes, "D_", model, "_0.1_qst_mtx.rds")))){
    gens <- data.matrix(read.table(paste0("input/fst0.1/geno_", demes, "D_", model, ".txt"), sep = ""))
    Qstmat <- readRDS(file.path(paste0("output/gaussian/", demes, "D_", model, "_0.1_qst_mtx.rds")))
    n <- rep(2000, demes)
    coal.mtx <- getMatrix(1000, demes, getPar(1000, demes, 0.1))
    Koch.mvn <- qst.sim(coal.mtx, nn=1e6)
    Nei.RoA <- RoA.Fst.Nei(gens, n)
    Nei.AoR <- AoR.Fst.Nei(gens, n)
    WC.RoA <- RoA.Fst.WC(gens, n)
    WC.AoR <- AoR.Fst.WC(gens, n)
    lknra.dens <- dchisq((demes - 1) * Fstat / Nei.RoA, df = demes - 1) * (demes - 1) / Nei.RoA
    lknar.dens <- dchisq((demes - 1) * Fstat / Nei.AoR, df = demes - 1) * (demes - 1) / Nei.AoR
    lkwra.dens <- dchisq((demes - 1) * Fstat / WC.RoA, df = demes - 1) * (demes - 1) / WC.RoA
    lkwar.dens <- dchisq((demes - 1) * Fstat / WC.AoR, df = demes - 1) * (demes - 1) / WC.AoR
    
    getLKMVN(1, model)
    getLKMVN(7, model)
    getLKMVN(2, model)
    getLKMVN(8, model)
  }
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", bty = "n", col = c('#88CCEE','#882255','#999933','#44AA99','#332288'), lwd = 2, horiz = TRUE, xpd = TRUE,
       legend = c('LK Nei RoA', 'LK Nei AoR', 'LK WC RoA', 'LK WC AoR', 'Koch MVN'),
       cex = 1.5)

dev.off()


for (demes in c(4,8,16)){
  print(paste0(demes, " demes"))
  pdf(paste0("figures/Figure_nulls_2_", demes, "D.pdf"),
      width = 10, height = 7.5)
  
  par(oma = c(2,1,1,1), mfrow=c(3,4), mar=c(6,6,1,1))
  for (model in c("split_star", "migration_island", "migration_circular")){
    print(model)
    if (any(file.exists(paste0("output/gaussian/", demes, "D_", model, "_0.1_qst_mtx.rds")))){
      gens <- data.matrix(read.table(paste0("input/fst0.1/geno_", demes, "D_", model, ".txt"), sep = ""))
      Qstmat <- readRDS(file.path(paste0("output/gaussian/", demes, "D_", model, "_0.1_qst_mtx.rds")))
      n <- rep(2000, demes)
      coal.mtx <- getMatrix(1000, demes, getPar(1000, demes, 0.1))
      Koch.mvn <- qst.sim(coal.mtx, nn=1e6)
      Nei.RoA <- RoA.Fst.Nei(gens, n)
      Nei.AoR <- AoR.Fst.Nei(gens, n)
      WC.RoA <- RoA.Fst.WC(gens, n)
      WC.AoR <- AoR.Fst.WC(gens, n)
      lknra.dens <- dchisq((demes - 1) * Fstat / Nei.RoA, df = demes - 1) * (demes - 1) / Nei.RoA
      lknar.dens <- dchisq((demes - 1) * Fstat / Nei.AoR, df = demes - 1) * (demes - 1) / Nei.AoR
      lkwra.dens <- dchisq((demes - 1) * Fstat / WC.RoA, df = demes - 1) * (demes - 1) / WC.RoA
      lkwar.dens <- dchisq((demes - 1) * Fstat / WC.AoR, df = demes - 1) * (demes - 1) / WC.AoR
      
      getLKMVN(1, model)
      getLKMVN(7, model)
      getLKMVN(2, model)
      getLKMVN(8, model)
    }
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend("bottom", bty = "n", col = c('#88CCEE','#882255','#999933','#44AA99','#332288'), lwd = 2, horiz = TRUE, xpd = TRUE,
         legend = c('LK Nei RoA', 'LK Nei AoR', 'LK WC RoA', 'LK WC AoR', 'Koch MVN'),
         cex = 1.5)
  
  dev.off()
}




######################################
# Supp Figure: Error rates
######################################
getER.all <- function(Qstvec){
  er1 <- colMeans(Qstvec > locusNei.cut)
  er2 <- colMeans(Qstvec > locusWC.cut)
  er3 <- colMeans(Qstvec > lknra.cut)
  er4 <- colMeans(Qstvec > lknar.cut)
  # er5 <- colMeans(Qstvec > lkwra.cut)
  # er6 <- colMeans(Qstvec > lkwar.cut)
  er7 <- colMeans(Qstvec > Mvn.cut)
  er8 <- colMeans(Qstvec > locusNeiCv.cut)
  # er9 <- colMeans(Qstvec > locusWCCv.cut)
  
  se1 <- sqrt(er1*(1 - er1)/nrow(Qstvec))
  lb1 <- er1 - 2*se1
  ub1 <- er1 + 2*se1
  se2 <- sqrt(er2*(1 - er2)/nrow(Qstvec))
  lb2 <- er2 - 2*se2
  ub2 <- er2 + 2*se2
  se3 <- sqrt(er3*(1 - er3)/nrow(Qstvec))
  lb3 <- er3 - 2*se3
  ub3 <- er3 + 2*se3
  se4 <- sqrt(er4*(1 - er4)/nrow(Qstvec))
  lb4 <- er4 - 2*se4
  ub4 <- er4 + 2*se4
  # se5 <- sqrt(er5*(1 - er5)/nrow(Qstvec))
  # lb5 <- er5 - 2*se5
  # ub5 <- er5 + 2*se5
  # se6 <- sqrt(er6*(1 - er6)/nrow(Qstvec))
  # lb6 <- er6 - 2*se6
  # ub6 <- er6 + 2*se6
  se7 <- sqrt(er7*(1 - er7)/nrow(Qstvec))
  lb7 <- er7 - 2*se7
  ub7 <- er7 + 2*se7
  se8 <- sqrt(er8*(1 - er8)/nrow(Qstvec))
  lb8 <- er8 - 2*se8
  ub8 <- er8 + 2*se8
  # se9 <- sqrt(er9*(1 - er9)/nrow(Qstvec))
  # lb9 <- er9 - 2*se9
  # ub9 <- er9 + 2*se9
  
  # ub.y <- max(er1,er2,er3,er4,er5,er6,er7,er8,er9)*5/4
  ub.y <- max(er1,er2,er3,er4,er7,er8)*5/4
  if (ub.y > 1){ub.y = 1}
  
  if (model == "split_star") {model.lb = "star-like"}
  if (model == "split_balanced") {model.lb = "balanced"}
  if (model == "split_caterpillar") {model.lb = "caterpillar"}
  if (model == "migration_island") {model.lb = "island"}
  if (model == "migration_circular") {model.lb = "circular"}
  
  par(mar=c(4,3,2,2))
  plot(x = log10(c(1, 10, 100, 1000)), y = c(0,0,0,0), ylim = c(0, ub.y), pch = "",
       bty = "n", ylab = "", xlab = "", xaxt = "n", yaxt = "n")
  axis(1, at = c(0, 1, 2, 3), labels = c("1", "10", "100", "1000"), cex.axis = 1.5)
  axis(2, las = 1, cex.axis = 1.5)
  points(log10(c(1, 10, 100, 1000)) - 0.05, er1[c(1,3,5,7)], col = pal[1], pch = 19)
  points(log10(c(1, 10, 100, 1000)) - 0.03, er2[c(1,3,5,7)], col = pal[2], pch = 19)
  points(log10(c(1, 10, 100, 1000)) - 0.01, er3[c(1,3,5,7)], col = pal[3], pch = 19)
  points(log10(c(1, 10, 100, 1000)) + 0.01, er4[c(1,3,5,7)], col = pal[4], pch = 19)
  # points(log10(c(1, 10, 100, 1000)) - 0.00, er5[c(1,3,5,7)], col = pal[5], pch = 19)
  # points(log10(c(1, 10, 100, 1000)) + 0.03, er6[c(1,3,5,7)], col = pal[6], pch = 19)
  points(log10(c(1, 10, 100, 1000)) + 0.03, er7[c(1,3,5,7)], col = pal[7], pch = 19)
  points(log10(c(1, 10, 100, 1000)) + 0.05, er8[c(1,3,5,7)], col = pal[8], pch = 19)
  # points(log10(c(1, 10, 100, 1000)) + 0.12, er9[c(1,3,5,7)], col = pal[9], pch = 19)
  points(log10(c(1, 10, 100, 1000)) - 0.05, er1[c(2,4,6,8)], col = pal[1], pch = 9)
  points(log10(c(1, 10, 100, 1000)) - 0.03, er2[c(2,4,6,8)], col = pal[2], pch = 9)
  points(log10(c(1, 10, 100, 1000)) - 0.01, er3[c(2,4,6,8)], col = pal[3], pch = 9)
  points(log10(c(1, 10, 100, 1000)) + 0.01, er4[c(2,4,6,8)], col = pal[4], pch = 9)
  # points(log10(c(1, 10, 100, 1000)) - 0.00, er5[c(2,4,6,8)], col = pal[5], pch = 9)
  # points(log10(c(1, 10, 100, 1000)) + 0.03, er6[c(2,4,6,8)], col = pal[6], pch = 9)
  points(log10(c(1, 10, 100, 1000)) + 0.03, er7[c(2,4,6,8)], col = pal[7], pch = 9)
  points(log10(c(1, 10, 100, 1000)) + 0.05, er8[c(2,4,6,8)], col = pal[8], pch = 9)
  # points(log10(c(1, 10, 100, 1000)) + 0.12, er9[c(2,4,6,8)], col = pal[9], pch = 9)
  lines(c(0,100), c(.05, .05), lty = 1.5, col = "grey")
  title(xlab = "Number of causal loci", line=2.5, cex.lab=1.4)
  title(bquote(.(demes)~D~.(model.lb)), adj = 1, line = 0, cex.lab=1.8)

  for(i in 1:8){
    lines(rep(floor((i - 1)/2) - .05, 2), c(lb1[i], ub1[i]), col = pal[1], lwd = 2)
    lines(rep(floor((i - 1)/2) - .03, 2), c(lb2[i], ub2[i]), col = pal[2], lwd = 2)
    lines(rep(floor((i - 1)/2) - .01, 2), c(lb3[i], ub3[i]), col = pal[3], lwd = 2)
    lines(rep(floor((i - 1)/2) + .01, 2), c(lb4[i], ub4[i]), col = pal[4], lwd = 2)
    # lines(rep(floor((i - 1)/2) - .00, 2), c(lb5[i], ub5[i]), col = pal[5], lwd = 2)
    # lines(rep(floor((i - 1)/2) + .03, 2), c(lb6[i], ub6[i]), col = pal[6], lwd = 2)
    lines(rep(floor((i - 1)/2) + .03, 2), c(lb7[i], ub7[i]), col = pal[7], lwd = 2)
    lines(rep(floor((i - 1)/2) + .05, 2), c(lb8[i], ub8[i]), col = pal[8], lwd = 2)
    # lines(rep(floor((i - 1)/2) + .12, 2), c(lb9[i], ub9[i]), col = pal[9], lwd = 2)
  }
}


pdf(paste0("figures/Figure_S14.pdf"), width = 8, height = 10)

par(oma = c(1,4,1,1), mfrow=c(4,3), mar=c(6,1,1,1))
for (demes in c(16,8,4,2)){
  print(paste0(demes, " demes"))
  for (model in c("split_star", "split_balanced", "split_caterpillar")){
    if (any(file.exists(paste0("output/gaussian/", demes, "D_", model, "_0.1_qst_mtx.rds")))){
      print(model)
      
      gens <- data.matrix(read.table(paste0("input/fst0.1/geno_", demes, "D_", model, ".txt"), sep = ""))
      Qstmat <- readRDS(file.path(paste0("output/gaussian/", demes, "D_", model, "_0.1_qst_mtx.rds")))
      n <- rep(2000, demes)
      Fstat <- seq(0,1, by = .0001)

      Fsts.Nei <- locus.Fsts.Nei(gens, n)
      Fsts.WC <- locus.Fsts.WC(gens, n)

      Nei.RoA <- RoA.Fst.Nei(gens, n)
      Nei.AoR <- AoR.Fst.Nei(gens, n)
      WC.RoA <- RoA.Fst.WC(gens, n)
      WC.AoR <- AoR.Fst.WC(gens, n)

      coal.mtx <- getMatrix(1000, demes, getPar(1000, demes, 0.1))
      Koch.mvn <- qst.sim(coal.mtx, nn=1e6)

      gens.freq <- cbind(gens, "freq" = sample.af.t(gens))
      gens.cv <- gens.freq[which(gens.freq[,ncol(gens.freq)]>0.05),-1]
      Fsts.Nei.cv <- locus.Fsts.Nei(gens.cv, n)
      Fsts.WC.cv <- locus.Fsts.WC(gens.cv, n)

      locusNei.cut <- quantile(Fsts.Nei, 0.95, na.rm = TRUE)
      locusWC.cut <- quantile(Fsts.WC, 0.95, na.rm = TRUE)
      lknra.cut <- qchisq(0.95, demes - 1) * Nei.RoA / (demes - 1)
      lknar.cut <- qchisq(0.95, demes - 1) * Nei.AoR / (demes - 1)
      lkwra.cut <- qchisq(0.95, demes - 1) * WC.RoA / (demes - 1)
      lkwar.cut <- qchisq(0.95, demes - 1) * WC.AoR / (demes - 1)
      Mvn.cut <- quantile(Koch.mvn, 0.95, na.rm = TRUE)
      locusNeiCv.cut <- quantile(Fsts.Nei.cv, 0.95, na.rm = TRUE)
      locusWCCv.cut <- quantile(Fsts.WC.cv, 0.95, na.rm = TRUE)
      
      getER.all(Qstmat)
    }
  }
}
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n', ann=FALSE)
legend("topleft", bty = "n", col = c(pal[c(1,2,3,4,7,8)], 'black', 'black'), pch = c(rep(19,7),9),
       legend = c(nulls[c(1,2,3,4,7,8)], expression("RB"~Q[ST]), expression("PBS"~Q[ST])), cex = 1.3)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
title(ylab="Type I error rate", line=-2, cex.lab=2.5)

dev.off()


pdf(paste0("figures/Figure_S15.pdf"), width = 10, height = 5)

par(oma = c(1,4,1,1), mfrow=c(2,4), mar=c(6,1,1,1))
for (model in c("migration_island", "migration_circular")){
  print(model)
  for (demes in c(16,8,4,2)){
    if (any(file.exists(paste0("output/gaussian/", demes, "D_", model, "_0.1_qst_mtx.rds")))){
      print(paste0(demes, " demes"))

      gens <- data.matrix(read.table(paste0("input/fst0.1/geno_", demes, "D_", model, ".txt"), sep = ""))
      Qstmat <- readRDS(file.path(paste0("output/gaussian/", demes, "D_", model, "_0.1_qst_mtx.rds")))
      n <- rep(2000, demes)
      Fstat <- seq(0,1, by = .0001)
      
      Fsts.Nei <- locus.Fsts.Nei(gens, n)
      Fsts.WC <- locus.Fsts.WC(gens, n)
      
      Nei.RoA <- RoA.Fst.Nei(gens, n)
      Nei.AoR <- AoR.Fst.Nei(gens, n)
      WC.RoA <- RoA.Fst.WC(gens, n)
      WC.AoR <- AoR.Fst.WC(gens, n)
      
      coal.mtx <- getMatrix(1000, demes, getPar(1000, demes, 0.1))
      Koch.mvn <- qst.sim(coal.mtx, nn=1e6)
      
      gens.freq <- cbind(gens, "freq" = sample.af.t(gens))
      gens.cv <- gens.freq[which(gens.freq[,ncol(gens.freq)]>0.05),-1]
      Fsts.Nei.cv <- locus.Fsts.Nei(gens.cv, n)
      Fsts.WC.cv <- locus.Fsts.WC(gens.cv, n)
      
      locusNei.cut <- quantile(Fsts.Nei, 0.95, na.rm = TRUE)
      locusWC.cut <- quantile(Fsts.WC, 0.95, na.rm = TRUE)
      lknra.cut <- qchisq(0.95, demes - 1) * Nei.RoA / (demes - 1)
      lknar.cut <- qchisq(0.95, demes - 1) * Nei.AoR / (demes - 1)
      lkwra.cut <- qchisq(0.95, demes - 1) * WC.RoA / (demes - 1)
      lkwar.cut <- qchisq(0.95, demes - 1) * WC.AoR / (demes - 1)
      Mvn.cut <- quantile(Koch.mvn, 0.95, na.rm = TRUE)
      locusNeiCv.cut <- quantile(Fsts.Nei.cv, 0.95, na.rm = TRUE)
      locusWCCv.cut <- quantile(Fsts.WC.cv, 0.95, na.rm = TRUE)
      
      getER.all(Qstmat)
    }
  }
}
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n', ann=FALSE)
legend("topleft", bty = "n", col = c(pal[c(1,2,3,4,7,8)], 'black', 'black'), pch = c(rep(19,7),9),
       legend = c(nulls[c(1,2,3,4,7,8)], expression("RB"~Q[ST]), expression("PBS"~Q[ST])), cex = 1.3)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
title(ylab="Type I error rate", line=-2, cex.lab=2.5)

dev.off()




######################################
# Supp Figure: Error rates across es
######################################
getER.es <- function(Qstvec){
  er1 <- colMeans(Qstvec > locusNei.cut)
  er2 <- colMeans(Qstvec > locusWC.cut)
  er3 <- colMeans(Qstvec > lknra.cut)
  er4 <- colMeans(Qstvec > lknar.cut)
  # er5 <- colMeans(Qstvec > lkwra.cut)
  # er6 <- colMeans(Qstvec > lkwar.cut)
  er7 <- colMeans(Qstvec > Mvn.cut)
  er8 <- colMeans(Qstvec > locusNeiCv.cut)
  # er9 <- colMeans(Qstvec > locusWCCv.cut)
  
  se1 <- sqrt(er1*(1 - er1)/nrow(Qstvec))
  lb1 <- er1 - 2*se1
  ub1 <- er1 + 2*se1
  se2 <- sqrt(er2*(1 - er2)/nrow(Qstvec))
  lb2 <- er2 - 2*se2
  ub2 <- er2 + 2*se2
  se3 <- sqrt(er3*(1 - er3)/nrow(Qstvec))
  lb3 <- er3 - 2*se3
  ub3 <- er3 + 2*se3
  se4 <- sqrt(er4*(1 - er4)/nrow(Qstvec))
  lb4 <- er4 - 2*se4
  ub4 <- er4 + 2*se4
  # se5 <- sqrt(er5*(1 - er5)/nrow(Qstvec))
  # lb5 <- er5 - 2*se5
  # ub5 <- er5 + 2*se5
  # se6 <- sqrt(er6*(1 - er6)/nrow(Qstvec))
  # lb6 <- er6 - 2*se6
  # ub6 <- er6 + 2*se6
  se7 <- sqrt(er7*(1 - er7)/nrow(Qstvec))
  lb7 <- er7 - 2*se7
  ub7 <- er7 + 2*se7
  se8 <- sqrt(er8*(1 - er8)/nrow(Qstvec))
  lb8 <- er8 - 2*se8
  ub8 <- er8 + 2*se8
  # se9 <- sqrt(er9*(1 - er9)/nrow(Qstvec))
  # lb9 <- er9 - 2*se9
  # ub9 <- er9 + 2*se9
  
  # ub.y <- max(er1,er2,er3,er4,er5,er6,er7,er8,er9)*5/4
  ub.y <- max(er1,er2,er3,er4,er7,er8)*5/4
  if (ub.y > 1){ub.y = 1}
  
  if (model == "split_star") {model.lb = "star-like"}
  if (model == "split_balanced") {model.lb = "balanced"}
  if (model == "split_caterpillar") {model.lb = "caterpillar"}
  if (model == "migration_island") {model.lb = "island"}
  if (model == "migration_circular") {model.lb = "circular"}
  
  par(mar=c(4,3,2,2))
  plot(x = log10(c(1, 10, 100, 1000)), y = c(0,0,0,0), ylim = c(0, ub.y), pch = "",
       bty = "n", ylab = "", xlab = "", xaxt = "n", yaxt = "n")
  axis(1, at = c(0, 1, 2, 3), labels = c("1", "10", "100", "1000"), cex.axis = 1.5)
  axis(2, las = 1, cex.axis = 1.5)
  points(log10(c(1, 10, 100, 1000)) - 0.05, er1[c(1,3,5,7)], col = pal[1], pch = 19)
  points(log10(c(1, 10, 100, 1000)) - 0.03, er2[c(1,3,5,7)], col = pal[2], pch = 19)
  points(log10(c(1, 10, 100, 1000)) - 0.01, er3[c(1,3,5,7)], col = pal[3], pch = 19)
  points(log10(c(1, 10, 100, 1000)) + 0.01, er4[c(1,3,5,7)], col = pal[4], pch = 19)
  # points(log10(c(1, 10, 100, 1000)) - 0.00, er5[c(1,3,5,7)], col = pal[5], pch = 19)
  # points(log10(c(1, 10, 100, 1000)) + 0.03, er6[c(1,3,5,7)], col = pal[6], pch = 19)
  points(log10(c(1, 10, 100, 1000)) + 0.03, er7[c(1,3,5,7)], col = pal[7], pch = 19)
  points(log10(c(1, 10, 100, 1000)) + 0.05, er8[c(1,3,5,7)], col = pal[8], pch = 19)
  # points(log10(c(1, 10, 100, 1000)) + 0.12, er9[c(1,3,5,7)], col = pal[9], pch = 19)
  points(log10(c(1, 10, 100, 1000)) - 0.05, er1[c(2,4,6,8)], col = pal[1], pch = 9)
  points(log10(c(1, 10, 100, 1000)) - 0.03, er2[c(2,4,6,8)], col = pal[2], pch = 9)
  points(log10(c(1, 10, 100, 1000)) - 0.01, er3[c(2,4,6,8)], col = pal[3], pch = 9)
  points(log10(c(1, 10, 100, 1000)) + 0.01, er4[c(2,4,6,8)], col = pal[4], pch = 9)
  # points(log10(c(1, 10, 100, 1000)) - 0.00, er5[c(2,4,6,8)], col = pal[5], pch = 9)
  # points(log10(c(1, 10, 100, 1000)) + 0.03, er6[c(2,4,6,8)], col = pal[6], pch = 9)
  points(log10(c(1, 10, 100, 1000)) + 0.03, er7[c(2,4,6,8)], col = pal[7], pch = 9)
  points(log10(c(1, 10, 100, 1000)) + 0.05, er8[c(2,4,6,8)], col = pal[8], pch = 9)
  # points(log10(c(1, 10, 100, 1000)) + 0.12, er9[c(2,4,6,8)], col = pal[9], pch = 9)
  lines(c(0,100), c(.05, .05), lty = 1.5, col = "grey")
  title(xlab = "Number of causal loci", line=2.5, cex.lab=1.4)
  title(bquote(.(es)~","~.(model.lb)), adj = 1, line = 0, cex.lab=1.8)
  
  for(i in 1:8){
    lines(rep(floor((i - 1)/2) - .05, 2), c(lb1[i], ub1[i]), col = pal[1], lwd = 2)
    lines(rep(floor((i - 1)/2) - .03, 2), c(lb2[i], ub2[i]), col = pal[2], lwd = 2)
    lines(rep(floor((i - 1)/2) - .01, 2), c(lb3[i], ub3[i]), col = pal[3], lwd = 2)
    lines(rep(floor((i - 1)/2) + .01, 2), c(lb4[i], ub4[i]), col = pal[4], lwd = 2)
    # lines(rep(floor((i - 1)/2) - .00, 2), c(lb5[i], ub5[i]), col = pal[5], lwd = 2)
    # lines(rep(floor((i - 1)/2) + .03, 2), c(lb6[i], ub6[i]), col = pal[6], lwd = 2)
    lines(rep(floor((i - 1)/2) + .03, 2), c(lb7[i], ub7[i]), col = pal[7], lwd = 2)
    lines(rep(floor((i - 1)/2) + .05, 2), c(lb8[i], ub8[i]), col = pal[8], lwd = 2)
    # lines(rep(floor((i - 1)/2) + .12, 2), c(lb9[i], ub9[i]), col = pal[9], lwd = 2)
  }
}


pdf(paste0("figures/Figure_S16.pdf"), width = 7.5, height = 7.5)

par(oma = c(1,4,1,1), mfrow=c(3,3), mar=c(6,1,1,1))
demes = 8
for (es in c("gaussian","uniform","laplace","alpha")){
  print(es)
  for (model in c("split_star","migration_island")){
    print(model)
    if (any(file.exists(paste0("output/gaussian/8D_", model, "_0.1_qst_mtx.rds")))){
      gens <- data.matrix(read.table(paste0("input/fst0.1/geno_8D_", model, ".txt"), sep = ""))
      Qstmat <- readRDS(file.path(paste0("output/", es, "/8D_", model, "_0.1_qst_mtx.rds")))
      n <- rep(2000, demes)
      
      Fsts.Nei <- locus.Fsts.Nei(gens, n)
      Fsts.WC <- locus.Fsts.WC(gens, n)
      
      Nei.RoA <- RoA.Fst.Nei(gens, n)
      Nei.AoR <- AoR.Fst.Nei(gens, n)
      WC.RoA <- RoA.Fst.WC(gens, n)
      WC.AoR <- AoR.Fst.WC(gens, n)
      
      coal.mtx <- getMatrix(1000, demes, getPar(1000, demes, 0.1))
      Koch.mvn <- qst.sim(coal.mtx, nn=1e6)
      
      gens.freq <- cbind(gens, "freq" = sample.af.t(gens))
      gens.cv <- gens.freq[which(gens.freq[,ncol(gens.freq)]>0.05),-1]
      Fsts.Nei.cv <- locus.Fsts.Nei(gens.cv, n)
      Fsts.WC.cv <- locus.Fsts.WC(gens.cv, n)
      
      locusNei.cut <- quantile(Fsts.Nei, 0.95, na.rm = TRUE)
      locusWC.cut <- quantile(Fsts.WC, 0.95, na.rm = TRUE)
      lknra.cut <- qchisq(0.95, demes - 1) * Nei.RoA / (demes - 1)
      lknar.cut <- qchisq(0.95, demes - 1) * Nei.AoR / (demes - 1)
      lkwra.cut <- qchisq(0.95, demes - 1) * WC.RoA / (demes - 1)
      lkwar.cut <- qchisq(0.95, demes - 1) * WC.AoR / (demes - 1)
      Mvn.cut <- quantile(Koch.mvn, 0.95, na.rm = TRUE)
      locusNeiCv.cut <- quantile(Fsts.Nei.cv, 0.95, na.rm = TRUE)
      locusWCCv.cut <- quantile(Fsts.WC.cv, 0.95, na.rm = TRUE)

      getER.es(Qstmat)
    }
  }
}
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n', ann=FALSE)
legend("topleft", bty = "n", col = c(pal[c(1,2,3,4,7,8)], 'black', 'black'), pch = c(rep(19,7),9),
       legend = c(nulls[c(1,2,3,4,7,8)], expression("RB"~Q[ST]), expression("PBS"~Q[ST])), cex = 1.3)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
title(ylab="Type I error rate", line=-2, cex.lab=2.5)

dev.off()




######################################
# Supp Table: Error rates for all
######################################
data.er <- matrix(NA, 0, 18)
for (model in c("split_star", "split_balanced", "split_caterpillar", "migration_island", "migration_circular")){
  print(model)
  for (demes in c(2,4,8,16)){
    if (any(file.exists(paste0("output/gaussian/", demes, "D_", model, "_0.1_qst_mtx.rds")))){
      print(paste0(demes, " demes"))
      
      gens <- data.matrix(read.table(paste0("input/fst0.1/geno_", demes, "D_", model, ".txt"), sep = ""))
      Qstmat <- readRDS(file.path(paste0("output/gaussian/", demes, "D_", model, "_0.1_qst_mtx.rds")))
      n <- rep(2000, demes)
      Fstat <- seq(0,1, by = .0001)
      
      Fsts.Nei <- locus.Fsts.Nei(gens, n)
      Fsts.WC <- locus.Fsts.WC(gens, n)
      
      Nei.RoA <- RoA.Fst.Nei(gens, n)
      Nei.AoR <- AoR.Fst.Nei(gens, n)
      WC.RoA <- RoA.Fst.WC(gens, n)
      WC.AoR <- AoR.Fst.WC(gens, n)
      
      coal.mtx <- getMatrix(1000, demes, getPar(1000, demes, 0.1))
      Koch.mvn <- qst.sim(coal.mtx, nn=1e6)
      
      gens.freq <- cbind(gens, "freq" = sample.af.t(gens))
      gens.cv <- gens.freq[which(gens.freq[,ncol(gens.freq)]>0.05),-1]
      Fsts.Nei.cv <- locus.Fsts.Nei(gens.cv, n)
      Fsts.WC.cv <- locus.Fsts.WC(gens.cv, n)
      
      locusNei.cut <- quantile(Fsts.Nei, 0.95, na.rm = TRUE)
      locusWC.cut <- quantile(Fsts.WC, 0.95, na.rm = TRUE)
      lknra.cut <- qchisq(0.95, demes - 1) * Nei.RoA / (demes - 1)
      lknar.cut <- qchisq(0.95, demes - 1) * Nei.AoR / (demes - 1)
      lkwra.cut <- qchisq(0.95, demes - 1) * WC.RoA / (demes - 1)
      lkwar.cut <- qchisq(0.95, demes - 1) * WC.AoR / (demes - 1)
      Mvn.cut <- quantile(Koch.mvn, 0.95, na.rm = TRUE)
      locusNeiCv.cut <- quantile(Fsts.Nei.cv, 0.95, na.rm = TRUE)
      locusWCCv.cut <- quantile(Fsts.WC.cv, 0.95, na.rm = TRUE)
      
      er1 <- colMeans(Qstmat > locusNei.cut)
      er2 <- colMeans(Qstmat > locusWC.cut)
      er3 <- colMeans(Qstmat > lknra.cut)
      er4 <- colMeans(Qstmat > lknar.cut)
      er5 <- colMeans(Qstmat > lkwra.cut)
      er6 <- colMeans(Qstmat > lkwar.cut)
      er7 <- colMeans(Qstmat > Mvn.cut)
      er8 <- colMeans(Qstmat > locusNeiCv.cut)
      er9 <- colMeans(Qstmat > locusWCCv.cut)
      
      data.er <- rbind(data.er, c(er1[7:8],er2[7:8],er3[7:8],er4[7:8],er5[7:8],er6[7:8],er7[7:8],er8[7:8],er9[7:8]))
    }
  }
}

data.er.rb <- data.er[,c(1,3,5,7,9,11,13,15,17)]
data.er.pbs <- data.er[,c(2,4,6,8,10,12,14,16,18)]
rownames(data.er.rb) <- c("2D star-like", "4D star-like", "8D star-like", "16D star-like",
                       "4D balanced", "8D balanced", "16D balanced",
                       "4D caterpillar", "8D caterpillar", "16D caterpillar",
                       "2D island", "4D island", "8D island", "16D island",
                       "4D circular", "8D circular", "16D circular")
colnames(data.er.rb) <- nulls
rownames(data.er.pbs) <- rownames(data.er.rb)
colnames(data.er.pbs) <- nulls

saveRDS(data.er.rb, file.path(paste0("figures/raw/Table_2_error_rates_QstRB.rds")))
write.csv(data.er.rb, "figures/raw/Table_2_error_rates_QstRB.csv")
saveRDS(data.er.pbs, file.path(paste0("figures/raw/Table_2_error_rates_QstPBS.rds")))
write.csv(data.er.pbs, "figures/raw/Table_2_error_rates_QstPBS.csv")


