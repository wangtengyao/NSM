#### Table 1 (Degree-based detection) ####
# this script generates all data in Table 1.

# load relevant funtions for simulation
source('algorithms.R')

n <- 500; N <- 30; tau <- 200; sd <- 1
ks <- c(3, 10, 30)
cs <- seq(0.05, 0.25, by=0.05)
nreps <- 100
old_mean <- 0.3
shape <- 10

#### fixed signed symmetric change (panel S1) ####
results <- array(0, c(length(ks), length(cs), nreps))

for (i in seq_along(ks)){
  for (j in seq_along(cs)){
    k <- ks[i]; c <- cs[j]
    for (rep in 1:nreps){
      set.seed(rep)
      x <- generateData(n, N, k, c, tau, shape=shape, symmetric=T, random_sign=F, old_mean=old_mean)$x_deg
      tau_hat <- findChange(x)$cp
      if (is.na(tau_hat)) {
        set.seed(rep)
        x <- generateData(n=2000, N, k, c, tau, shape=shape, symmetric=T, random_sign=F, old_mean=old_mean)$x_deg
        tau_hat <- findChange(x)$cp
      }
      results[i,j,rep] <- tau_hat - tau
      println('k=', k, ' c=', c, ' seed=', rep, ' declare=', tau_hat)
    }
  }
}

mean_false_alarm <- apply(setNA(results < 0, FALSE), c(1,2), mean)
mean_delay <- apply(results, c(1,2), function(x){mean(x[setNA(x >= 0, FALSE)])})
print(mean_delay)
df <- cbind(paste0('k=', ks), mean_false_alarm[,1], as.data.frame(mean_delay))
colnames(df) <- c(' ', 'false alarm rate', cs)
cat('\n\n\nTable 1 for Fixed Signed Symmetric Change:\n\n')
write.latextable(df, s=3)



#### fixed signed asymmetric change  (panel S2) ####
results <- array(0, c(length(ks), length(cs), nreps))

for (i in seq_along(ks)){
  for (j in seq_along(cs)){
    k <- ks[i]; c <- cs[j]
    for (rep in 1:nreps){
      set.seed(rep)
      x <- generateData(n, N, k, c, tau, shape=shape, symmetric=F, random_sign=F, old_mean=old_mean)$x_deg
      tau_hat <- findChange(x)$cp
      if (is.na(tau_hat)) {
        set.seed(rep)
        x <- generateData(n=2000, N, k, c, tau, shape=shape, symmetric=F, random_sign=F, old_mean=old_mean)$x_deg
        tau_hat <- findChange(x)$cp
      }
      results[i,j,rep] <- tau_hat - tau
      println('k=', k, ' c=', c, ' seed=', rep, ' declare=', tau_hat)
      #printPercentage(c(i,j,rep), dim(results))
    }
  }
}

mean_false_alarm <- apply(setNA(results < 0, FALSE), c(1,2), mean)
mean_delay <- apply(results, c(1,2), function(x){mean(x[setNA(x >= 0, FALSE)])})
print(mean_delay)
df <- cbind(paste0('k=', ks), mean_false_alarm[,1], as.data.frame(mean_delay))
colnames(df) <- c(' ', 'false alarm rate', cs)
cat('\n\n\nTable 1 for Fixed Signed Asymmetric Change:\n\n')
write.latextable(df, s=3)



#### random signed asymmetric change (panel S3) ####
results <- array(0, c(length(ks), length(cs), nreps))

for (i in seq_along(ks)){
  for (j in seq_along(cs)){
    k <- ks[i]; c <- cs[j]
    for (rep in 1:nreps){
      set.seed(rep)
      x <- generateData(n, N, k, c, tau, shape=shape, symmetric=F, random_sign=T, old_mean=old_mean)$x_deg
      tau_hat <- findChange(x)$cp
      if (is.na(tau_hat)) {
        set.seed(rep)
        x <- generateData(n=2000, N, k, c, tau, shape=shape, symmetric=F, random_sign=T, old_mean=old_mean)$x_deg
        tau_hat <- findChange(x)$cp
      }
      results[i,j,rep] <- tau_hat - tau
      println('k=', k, ' c=', c, ' seed=', rep, ' declare=', tau_hat)
      #printPercentage(c(i,j,rep), dim(results))
    }
  }
}

mean_false_alarm <- apply(setNA(results < 0, FALSE), c(1,2), mean)
mean_delay <- apply(results, c(1,2), function(x){mean(x[setNA(x >= 0, FALSE)])})
print(mean_delay)
df <- cbind(paste0('k=', ks), mean_false_alarm[,1], as.data.frame(mean_delay))
colnames(df) <- c(' ', 'false alarm rate', cs)
cat('\n\n\nTable 1 for Random Signed Asymmetric Change:\n\n')
write.latextable(df, s=3)
