# load packages
library(devtools)
library(putils)  # devtools::install_github('wangtengyao/putils')
library(InspectChangepoint)
library(igraph)

data=read.csv(file = 'Appended_symbol_similarity_30.csv')
weekIndex=read.csv('WeekIndex.csv')


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


N=30
n = round(nrow(data)/N,0)-3
Input=array(0, c(N,N,n))
x=matrix(0,N,n)
alpha = 0.01
lambda <- sqrt(2*log(n^2*N/alpha))
Date=weekIndex[1:(n-2),2]

for (i in 1:n)
{
  Input[1:N,1:N,i]=as.matrix(data[(((i-1)*N)+(i+1)):((i*N)+i), 1:N])
  x[,i]=log(rowSums(Input[1:N,1:N,i]))
}


vhat <- matrix(0, N, n)
for (t in 2:n){
  x.cusum <- cusum.transform(x[,1:t,drop=F])
  if (max(abs(x.cusum)) >= lambda)
    vhat[, t-1] <- sparse.svd(x.cusum, lambda)
  else
    vhat[, t-1] = PiW(rnorm(N))  # random draw from Unit Sphere
  # printPercentage(t,2:n)
}

A = rep(0,n-2)
sinA <- rep(0, n-2)
for (t in 1:(n-2)){
  A[t] <- sinTheta(vhat[,t], vhat[,t+1])$angle
  sinA[t] <- sinTheta(vhat[,t], vhat[,t+1])$sinL
}

alpha=0.12
b=round(log(n)/(N/8-log(2/alpha)))
print(b)
t=b+1
while(mean(sinA[t:(t-b)])>0.5)
# while(sinA[t]>0.5 )
{ t=t+1}
print(t)
# [1] 30
vhat[, t]
binary_vhat=vhat
binary_vhat[which(binary_vhat!=0)]=1
binary_vhat[which(binary_vhat==0)]=NA

barplot(binary_vhat, main = "Coordinates govern projection direction", ylab = "number of coordinate", xlab = "Time(week)")
abline(v=35,col='green', lwd=3)

# plot(binary_vhat[,27],ylab=names(vhat[27]),type="p")
# for ( i in seq(1,ncol( binary_vhat ),1) ) plot(binary_vhat[,i],ylab=names(vhat[i]),type="p")


par(mfrow=c(1,2),mar=c(2,2,2,2))
A_week=data.frame(cbind(Date,A), index=Date)
plot(A_week$A~A_week$index,type="l",xlab="Date", main="angle")
sinA_week=data.frame(cbind(Date,sinA), index=Date)
plot(sinA_week$sinA~sinA_week$index,type="l",xlab="Date",main="sin(angle)")


########################################################################
######################## Graph/network Visulation  #####################
########################################################################


G=Input[,,5]
G=as.matrix(G)
diag(G)=0
G[which(G!=0)]=1
g5 <- graph.adjacency(G, mode="undirected")

G=Input[,,20]
G=as.matrix(G)
diag(G)=0
G[which(G!=0)]=1
g20 <- graph.adjacency(G, mode="undirected")


G=Input[,,30]
G=as.matrix(G)
diag(G)=0
G[which(G!=0)]=1
g29 <- graph.adjacency(G, mode="undirected")

G=Input[,,45]
G=as.matrix(G)
diag(G)=0
G[which(G!=0)]=1
g45 <- graph.adjacency(G, mode="undirected")


# png("network.png",width=2400,height=1600)
par(mfrow=c(2,2),mar=c(2,2,2,2))
# pal2 <- rainbow(5, alpha=.5)
plot(g5, pch=19, cex=5, col='blue', vertex.color="gold", vertex.size=20, main=Date[5])
plot(g20, pch=19, cex=5, col='blue', vertex.color="gold", vertex.size=20, main=Date[20])
plot(g30, pch=19, cex=5, col='blue', vertex.color="gold", vertex.size=20, main=Date[30])
plot(g45, pch=19, cex=5, col='blue', vertex.color="gold", vertex.size=20, main=Date[45])

# dev.off()






