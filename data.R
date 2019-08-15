#################standard normal example###################

numofsims=50;
n=300
p=2
modelchoice=2

data3=array(0,c(n,p,numofsims))
clust.true=matrix(0,numofsims,n)
for (simnum in 1:numofsims){
  data3[,,simnum]=rmvnorm(n,c(0,0),diag(1,2,2))
  clust.true[simnum,]=matrix(1,1,n)
}


df = data.frame(x=data3,cluster=clust.true)
ggplot(data=df, aes(x=df$x.1, y=df$x.2)) + 
  geom_point() +
  ggtitle("Normal Data") +xlab("x1") +ylab("x2")



################uniform example##########################
numofsims=10;

n=100
p=2

data3=array(0,c(n,p,numofsims))
clust.true=matrix(0,numofsims,n)
for (simnum in 1:numofsims){
  rsim=sqrt(runif(n))
  thetasim=runif(n)*2*pi
  data3[,,simnum]=cbind(rsim*sin(thetasim),rsim*cos(thetasim))
  clust.true[simnum,]=matrix(1,1,n)
}

# Plot simulation
simnum=1
par(mfrow=c(1,1),mar=c(2.5,2.5,.5,.5)+.1, mgp=c(1.5, 0.5, 0))
plot(data3[,1,simnum],data3[,2,simnum],xlab="x_1",ylab="x_2",main="Uniform",xlim=c(-1.1,1.1),ylim=c(-1.1,1.1))
xgrid=seq(-1,1,.01) 
lines(xgrid,sqrt(1-xgrid^2))
lines(xgrid,-sqrt(1-xgrid^2))




#################mixture 2 bivariate normal example###################
library(ggplot2)
library(mvtnorm)



numofsims=50;

n=100
p=2
modelchoice=2
data3=array(0,c(n,p,numofsims))
clust.true=matrix(0,numofsims,n)
for (simnum in 1:numofsims){
  aux=runif(n)
  clust.true[simnum,aux<0.3]=1
  clust.true[simnum,aux>0.3]=2
  data3[clust.true[simnum,]==1,,simnum]=rmvnorm(sum(clust.true[simnum,]==1),c(0,0),diag(1,2,2))
  data3[clust.true[simnum,]==2,,simnum]=rmvnorm(sum(clust.true[simnum,]==2),c(3,3),diag(1,2,2))

}


df <- data.frame(mu=data3[,,1],cluster=clust.true[1,])

ggplot(df, aes(x=df$mu.1,y=df$mu.2)) + geom_point(aes(group=cluster, colour=cluster))



#################mixture 3 univariate normal example###################

library(ggplot2)

n=100
aux=runif(n)
clust.true=rep(0,n)
clust.true[aux<0.3]=1
clust.true[aux>0.3& aux<0.6]=2
clust.true[aux>0.6]=3
data4=matrix(0,n,1)
data4[clust.true==1,]=rnorm(sum(clust.true==1),0,1)
data4[clust.true==2,]=rnorm(sum(clust.true==2),5,1)
data4[clust.true==3,]=rnorm(sum(clust.true==3),10,1)

