setwd("")

rm(list = ls())

library(diversitree)

## Read the phylogenetic trees

tr <- read.nexus("tree.nex")

trs <- read.tree("trees.nex")

tr100 <- sample(trs, 100)

# function to fix subplex error
fixerror <- function(lik2, p) {
  res <- try(find.mle(lik2, p[argnames(lik2)], method = "subplex"))
  if ( inherits(res, "try-error") ) 
    res <- list(lnLik = -Inf)
  res
}

##### Step one #####

library(phytools)

## Simulate a Character Distribution on a Tree

# Fit a trait-independent model to the data (e.g., the mkn model, in this case)
# Read character states data, is the same matrix of the GeoSSE analysis
char <- read.csv("states.csv", row.names = 1)

# Character data for analysis
sp <- char[trs[[1]]$tip.label, ]
names(sp) <- trs[[1]]$tip.label

#Using ACE
fitER <- ace(sp, tr, model = "ER", type = "discrete")
fitER$index.matrix

round(fitER$lik.anc, 3)

#Using SIMMAP
mtree <- make.simmap(tr, sp, model = "ER")

q <- mtree$Q
q
Q <- q
colnames(Q) <- c(1, 0, 2)
rownames(Q) <- c(1, 0, 2)
Q

# Simulate the traits using the matrix of transitions of SIMMAP
sim <- phytools::sim.history(tr, Q, anc = NULL, nsim = 1)$states #trait simulation

simutrait <- list()

for(i in 1:100){
  print(i)
  simutrait[[i]] <- sim.history(tr100[[i]], Q, anc = NULL)$states
  #matfur <- data.frame(matrix(unlist(resu), nrow = 652, byrow=T))
  #traitsim <- cbind(as.matrix(taxnames), matfur)
}

simutrait[[1]]

# Using lapply, this return the list object of simulated character states

z <- as.numeric(traitsim[[1]])
names(z) <- tr$tip.label
z

foo <- function(tree, Q) sim.history(tree, Q, anc = NULL)$states

traitsim <- lapply(tr100, foo, Q)

sims <- lapply(traitsim, as.numeric)

simsl <- list() #assign species names to the list of data
for (i in 1:100){
  print(i)
   foo <- sims[[i]]
   names(foo) <- tr$tip.label
   simsl[[i]] <- foo 
   }

simsl[[2]]

##### Step 2 #####

library(diversitree)

# Simulation based GeoSSE
sim.full.lik <- list()

sim.null.lik <- list()

sim.fit.full <- list()

sim.fit.null <- list()

full.lik <- list()
null.lik <- list()
full.fit <- list()
null.fit <- list()

for(i in 1:100){
  
  print(i)
  
  ltrs <- trs[[i]]

sim.full.lik <- lapply(simsl, make.geosse, tree = ltrs) # 'lapply' just
# applies 'make.geosse' to all the elements of 'sim.dis'
p <- starting.point.geosse(ltrs)

full.lik[[i]] <- sim.full.lik

sim.null.lik <- lapply(sim.full.lik, constrain, 
                       formulae = list("sA ~ sB", "xA ~ xB", "dA ~ dB", "sAB ~ 0"))
null.lik[[i]] <- sim.null.lik
}

sim.fit.full <- lapply(unlist(full.lik), find.mle, p, method = 'subplex')
sim.fit.full
sim.fit.null <- lapply(unlist(null.lik), find.mle, p[-c(2, 3, 5, 7)])
sim.fit.null

null.dist <- numeric(length = 10000)

for (i in 1:10000) null.dist[i] <- 2 * (sim.fit.full[[i]]$lnLik -
                                          sim.fit.null[[i]]$lnLik)
hist(null.dist)# See what the null distribution looks like
lr.emp <- 2 * (fit.full$lnLik - fit.null$lnLik)
#lr.emp <- 2 * (-2298.3115 - -2363.5071)
abline(v = lr.emp, col = 'red') # Show where the empirical LR falls
mean(lr.emp <= null.dist) # Simulation-based P-value

save.image("SimulationFinal.RData")

##### Observed and simulated  models #####

full <- read.csv("Observed_best.csv")
null <- read.csv("Observed_null.csv")

dLL <- 2 * (full$lnLik - null$lnLik)

meandLL <- mean(dLL)

hist(dLL)
abline(v = meandLL, col = "red")

hist(dLL, col = "gray", add = T, lty = 0)

##### Final plot of the simulation-based p-value with an inset of the observed dLL
library(Hmisc)

tiff(filename = "p_value_inset_final.tif", width = 15, height = 20, units = "cm", 
     pointsize = 12, compression = "lzw", bg = "white", res = 600)

par(oma = c(2, 2.5, 1, 1))
par(mar = c(3, 4, 1, 1))
par(lwd = 1.5)

hist(null.dist, col = "black", main = NULL, xlim = c(0, 300), cex.axis = 3,
     xlab = NULL, ylab = NULL, lty = 0, yaxt = "n", xaxt = "n"); axis(2, labels = F, 
                                                                      cex.axis = 1.3)
axis(1, labels = T, cex.axis = 1.3)
abline(v = meandLL, col = "darkgray", lwd = 3, lty = 2)
title(xlab = "Simulated ??LL", outer = T, line = 0, cex = 3, font = 2, cex.lab = 1.5)
title(ylab = "Frequency", outer = T, line = -1, cex = 3, font = 2, cex.lab = 1.5)
box()
inset <- subplot(hist(dLL, main = NULL, xlab = NULL, ylab = NULL, 
                      yaxt = "n", xlim = c(50, 300)), 
        x = 243, y = 1177.5, size = c(1.8, 1.8)); axis(2, labels = F)

op <- par(no.readonly = TRUE)
par(inset)
abline(v = meandLL, col = "darkgray", lwd = 3, lty = 2)
par(op)
op1 <- par(no.readonly = TRUE)
par(inset)
box()
par(op1)

dev.off()

