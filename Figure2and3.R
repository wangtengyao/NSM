##### This script generates all figures in Fig 2 and 3 #####

# load relevant simulation functions
source('algorithms.R')

###############################################################
##############  change point test for edges  #################
###############################################################


# model parameters
N <- 30
n <- 500
tau <- 200
old_mean <- 0.3
shape <- 10

# sinTheta function

sinTheta <- function(U,V)
{
  output=NULL
  angle=acos((t(U)%*%V)/(vector.norm(U)*vector.norm(V)))
  output$angle=angle
  sinL=sin(angle)
  output$sinL=sinL
  return(output)
}


########################################################################
######################## Visulation for simulation #####################
########################################################################

visualiseMethod <- function(k, c, symmetric, random_sign, figure_output){
  # generate data for change point simulation
  x <- generateData(n, N, k, c, tau, shape, symmetric, random_sign, old_mean)$x_deg
  # compute the sequence of projection vectors in vhat and detect change using
  # the findChange() function
  tmp <- findChange(x)
  cp <- tmp$cp; vhat <- tmp$vhat

  # compute angle between consecutive projections and the sine of the angles
  A <- rep(NA, n)
  for (t in 3:n){
    A[t] <- acos(sum(vhat[,t-1] * vhat[,t]))
  }
  sinA <- sin(A)

  # plot the figure
  png(figure_output, width=800, height=600)
  palet <- matplotlib_palette(10)
  par(mfrow=c(1,2),mar=c(2,2,2,2))
  plot(A, xlab='Date', ylab='angle', main='angle', col=palet[1])
  abline(v=c(tau, cp), col=palet[3:2], lty=1)

  plot(sinA, xlab='Date', ylab='angle', main='angle', col=palet[1])
  abline(v=c(tau, cp), col=palet[3:2], lty=1)
  dev.off()
}


##### main (Figure 2) #####

# 2a: sparse change, small signal
visualiseMethod(k=3, c=0.05, symmetric=TRUE, random_sign=FALSE,
                figure_output='Fig2a.png')

# 2b: dense change, small signal
visualiseMethod(k=30, c=0.05, symmetric=TRUE, random_sign=FALSE,
                figure_output='Fig2b.png')

# 2c: sparse change, big signal
visualiseMethod(k=3, c=0.25, symmetric=TRUE, random_sign=FALSE,
                figure_output='Fig2c.png')

# 2d: dense change, big signal
visualiseMethod(k=30, c=0.25, symmetric=TRUE, random_sign=FALSE,
                figure_output='Fig2d.png')




##### main (Figure 3) #####

# 3a: sparse change, small signal
visualiseMethod(k=3, c=0.05, symmetric=FALSE, random_sign=TRUE,
                figure_output='Fig3a.png')

# 3b: dense change, small signal
visualiseMethod(k=30, c=0.05, symmetric=FALSE, random_sign=TRUE,
                figure_output='Fig3b.png')

# 3c: sparse change, big signal
visualiseMethod(k=3, c=0.25, symmetric=FALSE, random_sign=TRUE,
                figure_output='Fig3c.png')

# 3d: dense change, big signal
visualiseMethod(k=30, c=0.25, symmetric=FALSE, random_sign=TRUE,
                figure_output='Fig3d.png')



