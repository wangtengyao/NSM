source('algorithms.R')
library("foreach")
library("doParallel")
library("InspectChangepoint")
library("putils")

n <- 500; N <- 30; tau <- 200;
c <- 0.02

ks <- c(2,3,5) # changed nodes
ksE <- c(10,15,20,25,30)   # changed edges per node
#cs <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5)
cs <- c(0.02,0.05,0.1, 0.15, 0.2, 0.5)


################################
################################
########### edges ##############
################################
################################

nreps <- 128
resultsEdges <- array(0, c(length(ks), length(ksE), length(cs), 2, nreps))

generateDataSomeEdges <- function(n, N, kE, k, c, tau, shape=10, symmetric=TRUE, random_sign=TRUE){
  # compute mean edge weight matrix before and after change
  A0 <- matrix(0.3, N, N); diag(A0) <- 0 # mean edge weight matrix before change
  A1 <- A0

  # we make changes in the first kE edges of the first k nodes (k=2 so far, but we can let it vary too)
  change_edges <- outer(1:N, 1:N, function(i,j){i<=k & i<j & j<=kE})
  #change_edges <- outer(1:N, 1:N, function(i,j){i<=k & i!=j & j<=kE})

  A1[change_edges] <- A0[change_edges] +
    c * sample(c(1,(-1)^random_sign), sum(change_edges), replace=TRUE)
  if (symmetric) A1[lower.tri(A1)] <- t(A1)[lower.tri(A1)]

  # generate edge weights and compute degree vectors
  if (symmetric) {xE <- matrix(0, N*(N-1)/2, n)} else {xE <- matrix(0, N^2-N, n)}
  if (symmetric) {x <- matrix(0, N, n)} else {x <- matrix(0, 2*N, n)}

  for (t in 1:n){
    # generate edge weights from Beta distribution
    if (t<=tau) {mu <- A0} else {mu <- A1}
    E <- matrix(rbeta(N^2, shape*mu, shape*(1-mu)), N, N)

    if (symmetric) {
      E[lower.tri(E)] <- t(E)[lower.tri(E)]
      ##### YO:    we monitor both using nodes and edges; xE for edges and x for nodes
      xE[,t] <- E[lower.tri(E)];
      x[,t] <- rowSums(E)
    } else {
      xE[,t] <- c(E[lower.tri(E)], E[upper.tri(E)]) # out-degree followed by in-degree
      x[,t] <- c(rowSums(E), colSums(E)) # out-degree followed by in-degree
    }
  }
  return(list(x=x, xE=xE)); # we return both a sequence of edges and a sequence of degrees
}


findChange <- function(x, alpha=0.01, b=2, lambda_scale=1){
  N <- nrow(x); n <- ncol(x)
  # compute matrix of projection directions
  vhat <- matrix(0, N, n)
  for (t in 2:n){
    cat(t,"\n")
    lambda <- sqrt(2*log(t^2*N/alpha)) * lambda_scale
    y <- x[, 1:t, drop=F]
    y <- y / mad(y[,-1]-y[,-t]) * sqrt(2)  # normalise x by current estimate of sd
    x.cusum <- cusum.transform(y)
    #cat(max(abs(x.cusum)), "\n")
    if (max(abs(x.cusum)) >= lambda)  {vhat[, t] <- sparse.svd(x.cusum, lambda)}  else  {vhat[, t] <- random.UnitVector(N)}
  }

  # compute consecutive sine theta distance of projections
  A <- rep(1, n)
  for (t in 3:n){
    A[t] <- sinThetaLoss(vhat[, t-1], vhat[, t])
  }

  tmp <- A < 0.5
  for (t in b:n){
    if (all(tmp[(t-b+1):t])) return(t)
  }
  return(NA)
}

### YO: parallelized locally using foreach (repetions only)
cl <- makeCluster(detectCores())
cat("Cores", detectCores())
registerDoParallel(cl)



for (ik in 1:length(ks))
{
for (i in seq_along(ksE)){   # we change a subset of edges per node
  for (j in seq_along(cs)){   # size of the shift as before
    kE <- ksE[i]; c <- cs[j]; k <- ks[ik];
    cat("ik ", ik, " i ", i, " j ", j, "\n")

    results.temp <-   foreach(rep=1:nreps, .combine="cbind", .inorder=FALSE, .packages=c("InspectChangepoint", "putils")) %dopar%
    {
      set.seed(rep)

      x <- generateDataSomeEdges(n, N, kE, k, c, tau)
      tau_hat <- findChange(x$x); # find change using degrees
      tau_hatE <- findChange(x$xE); # find change using edges
      if (is.na(tau_hat)) {
        x <- generateDataSomeEdges(n=2000, N,kE, k, c, tau)
        tau_hat <- findChange(x$x)
      }
      if (is.na(tau_hatE)) {
        x <- generateDataSomeEdges(n=2000, N,kE, k, c, tau)
        tau_hatE <- findChange(x$xE)
      }
      cat(rep, "  ", tau_hat-tau, "  ", tau_hatE-tau, "\n");
      return(c(tau_hat-tau, tau_hatE-tau));   # return delay for the two monitoring approaches
    }

    resultsEdges[ik,i,j,,] <- results.temp;

      #results[i,j,rep] <- tau_hat - tau
      #printPercentage(c(i,j,rep), dim(results))
  }
}
}


save(resultsEdges, file="resultsSomeEdges_symmetric_k2-3-5.RData")
stopCluster(cl);

load(file="resultsSomeEdges_symmetric_k2-3-5.RData")

for (ik in 1:length(ks))
{

png(paste("edge_initiated_log_k", ks[ik], ".png", sep=""), width=600,height=500)
palette=putils::matplotlib_palette(10, visualise=T)
par(mar = c(5, 5, 2, 2))
matplot(cs[-1], t(mean_delay[ik,,-1,1]), type ="b", ylab = "Delay", log="y",
        xlab = "c", col = palette[1],lwd=2, pch=20, ylim=c(10,500))
matplot(cs[-1], t(mean_delay[ik,,-1,2]), type ="b", col = palette[2],lwd=2, pch=20, add=TRUE)
legend("topright", c("node detection", "edge detection"),   col = palette[1:2], lty=1, lwd=2, bty="n")
legend("bottomleft", legend=ksE,   col = palette[1], lty=1:5, lwd=2, bty="n")

dev.off()

}




mean_false_alarm <- apply(setNA(resultsEdges < 0, FALSE), c(1,2,3), mean)

mean_delay <- apply(resultsEdges, c(1,2,3,4), function(x){mean(x[setNA(x >= 0, FALSE)])})
print(mean_delay)
df <- cbind(paste0(ks/N*100, '\\%'), mean_false_alarm[,1], as.data.frame(mean_delay))
colnames(df) <- c(' ', 'false alarm rate', cs)
write.latextable(df, file='outputSomeEdges.txt', s=3)



