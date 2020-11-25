##### This R script generates figure 1 of the paper #####
##### it investigates the effect of the tuning parameter lambda on both
##### the false alarm rate and response delay in various simulation settings


# load relevant simulation functions #
source('algorithms.R')

##### main simulation function #####

FAvRD <- function(n=1000, N=30, tau=200, sd=1, nreps=100, old_mean=0.3,
                  k=10, shape=10, c=0.02, symmetric=T, random_sign=F)
{
  param_str <- paste0('k=',k,'_c=',c,'_sym=',symmetric,'_rand=',random_sign)
  lambda_scales <- exp(seq(log(0.5),log(1.5),length=20)) * sqrt(2)
  results <- matrix(0, length(lambda_scales), nreps)

  for (rep in 1:nreps){
    set.seed(rep)
    x <- generateData(n, N, k, c, tau, shape, symmetric, random_sign, old_mean)$x_deg
    for (i in seq_along(lambda_scales)){
      lambda_scale <- lambda_scales[i]
      results[i, rep] <- findChange(x, lambda_scale=lambda_scale/sqrt(2))$cp
      println(param_str, ': ', rep, ' ', i, ' ', results[i, rep])
    }
  }

  FA <- rowMeans(results <= tau)
  tmp <- results; tmp[tmp <= tau] <- NA
  RD <- rowMeans(tmp, na.rm=T) - tau
  println()
  println(param_str)
  print(rbind(FA=FA, RD=RD))
  println()

  save(results, file=paste0('FAvRD_', param_str, '.RData'))
}

##### function for generating figures #####
plotFigure <- function(data_input, figure_output){
  # load saved simulation data
  load(data_input)

  # compute FA and RD values to plot
  lambda_scales <- exp(seq(log(0.5),log(1.5),length=20)) * sqrt(2)
  tau <- 200
  FA <- rowMeans(results <= tau)
  tmp <- results; tmp[tmp <= tau] <- NA
  RD <- rowMeans(tmp, na.rm=T) - tau

  # create figure
  png(figure_output, width=500,height=350)
  palette=putils::matplotlib_palette(10)
  par(mar = c(5, 5, 3, 5))

  # plot False Alarm rate against scale parameter r in the definition of lambda
  plot(lambda_scales, FA, type ="b", ylab = "False alarm rate",
       xlab = "r", col = palette[1],lwd=2, pch=20)

  # plot Response Delay against the r sequence
  par(new = TRUE)
  plot(lambda_scales, RD, type = "b", xaxt = "n", yaxt = "n",
       ylab = "", xlab = "", ylim=c(0,40), col = palette[2], lwd=2, pch=20)
  axis(side = 4)
  mtext("Response delay", side = 4, line = 3)

  # draw our theoretical and practical suggestion for r value
  abline(v=c(sqrt(2),2), lty=2, col=palette[c(3,4)])

  # figure legend
  legend("topright", c("false alarm rate", "response delay"),
         col = palette[1:2], lty=1, lwd=2)
  dev.off()
}


##### main #####

# run simulations (can be skipped and use the uploaded .RData file directly)
FAvRD(k=3, c=0.05, symmetric=TRUE, random_sign=FALSE)
FAvRD(k=30, c=0.25, symmetric=FALSE, random_sign=TRUE)

# plot figures
plotFigure(data_input='FAvRD_k=3_c=0.05_sym=TRUE_rand=FALSE.RData',
           figure_output='Fig1a.png')
plotFigure(data_input='FAvRD_k=30_c=0.25_sym=FALSE_rand=TRUE.RData',
           figure_output='Fig1b.png')
