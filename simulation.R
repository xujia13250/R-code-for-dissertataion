setwd("C:/Users/Administrator/Desktop/ÂÛÎÄ/ÂÛÎÄ/part3/")
#install.packages("devtools")
#devtools::install_github('sarawade/mcclust.ext')
source("credibleball_v2.R")
source("gibbs3.R")
library(mvtnorm) 
library(mcclust)
library(msm)
library(MCMCpack)
library(mcclust.ext)
library(HDInterval)
library(ggplot2)
library(dirichletprocess)
############simulate the data(GMM)########################


## Initialise output
modelchoice=2
kVIsims=matrix(0,numofsims,modelchoice)
kBindersims=matrix(0,numofsims,modelchoice)
kMAPsims=matrix(0,numofsims,modelchoice)
VIsims=matrix(0,numofsims,modelchoice)
MAPsims=matrix(0,numofsims,modelchoice)
Bindersims=matrix(0,numofsims,modelchoice)
VIclus=array(0,c(numofsims,n,modelchoice))
MAPclus=array(0,c(numofsims,n,modelchoice))
Binderclus=array(0,c(numofsims,n,modelchoice))
khatsims=matrix(0,numofsims,modelchoice)
kmodesims=matrix(0,numofsims,modelchoice)
maphpd <- matrix(0,numofsims,modelchoice)
alphacb=c(0.01,seq(0.025,0.975,0.025),0.99)
cbVIeps=array(0,c(numofsims,length(alphacb),modelchoice))
cbBindereps=array(0,c(numofsims,length(alphacb),modelchoice))
cbMAPeps=array(0,c(numofsims,length(alphacb),modelchoice))

#####################################################################



################function########################
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


#DPMM.post <- function(x,result,niter,burnin){
  #data <- scale(x)
  #z <- plyr::ldply(result$labelsChain, rbind)[(burnin+1):(burnin+niter),]
  #loglik <- result$likelihoodChain[(burnin+1):(burnin+niter)]
  #alpha <- result$alphaChain[(burnin+1):(burnin+niter)]
  #alpha.para <- result$alphaPriorParameters
  #n <- dim(z)[1]
  #d <- dim(x)[2]
  #post.mat <- numeric(n)
  #for (i in 1:n){
  #  post.log <- numeric(1)
  # df <- data.frame(x=data,(z[i,]),row.names = NULL)
    #colnames(df)[3] <-'z'
    #mu <- result$clusterParametersChain[[burnin+i]]$mu
    #sigma <- result$clusterParametersChain[[burnin+i]]$sig
    #max.index <- length(unique(df$z))
    #for (j in 1:max.index){

   #   post.log <- post.log + sum(dmvnorm(mu[,,j],c(0,0),2*sigma[,,j],log=TRUE))
  #                       + log(dwish(sigma[,,j],2,diag(1,2,2))) 
                           
 #     }
#    post.mat[i] <- post.log + loglik[i]+ log(dgamma(alpha[burnin+i],alpha.para[1],alpha.para[2]))
#  }
 # return(post.mat)
#}

DPMM.post <- function(x,result,T0=diag(1,2,2),kappa0=2,v0=2,mu0=c(0,0),niter,burnin){
  data <- scale(x)
  z <- plyr::ldply(result$labelsChain, rbind)[(burnin+1):(burnin+niter),]
alpha <- result$alphaChain[(burnin+1):(burnin+niter)]
alpha.para <- result$alphaPriorParameters
  n <- dim(z)[1]
  d <- dim(x)[2]
  post.mat <- matrix(0,ncol=n)
  for (i in 1:n){
    post.log <- numeric(1)
    df <- data.frame(x=data,t(as.matrix(z[i,])),row.names = NULL)
    colnames(df)[3] <-'z'
    mu <- result$clusterParametersChain[[burnin+i]]$mu
    sigma <- result$clusterParametersChain[[burnin+i]]$sig
    max.index <- length(unique(df$z))

    for (j in 1:max.index){
      T_k <- matrix(0,d,d)
      y.diff <- sweep(df[df$z==j,][,1:d],2,colMeans(df[df$z==j,][,1:d]),'-')
      y.diff.sum <- matrix(0,d,d)
      for (h in 1:dim(y.diff)[1]){
        y.diff.sum <- y.diff.sum + t(y.diff[h,])%*%as.matrix(y.diff[h,])
      }
      T_k <- T0 + y.diff.sum + kappa0*length(z[i,]==j)/
        (kappa0+length(z[i,]==j))*
        as.matrix((mu0-colMeans(df[df$z==j,][,1:d])))%*%t((mu0-colMeans(df[df$z==j,][,1:d])))
      v=v0+length(z[i,]==j)
      mu_k <- (kappa0*mu0+length(z[i,]==j)*colMeans(df[df$z==j,][,1:d]))/(kappa0+length(z[i,]==j))
      kappa = kappa0 + length(z[i,]==j)
      post.log <- post.log + sum(dmvnorm(df[df$z==j,][,1:d],mean=mu[,,j], sigma=sigma[,,j],log = TRUE))+
      sum(dmvnorm(mu[,,j],mu_k,1/kappa*sigma[,,j],log=TRUE)) + log(diwish(sigma[,,j],v,T_k))
    }
    post.mat[1,i] <- post.log + log(dgamma(alpha[i],alpha.para[1],alpha.para[2]))
  }  
  return(post.mat)
}



gibbs.lpost.clusttrue = function(x,clust,T0=diag(1,2,2),kappa0=2,v0=2,mu0=c(0,0),d=2){
  
  df = data.frame(x=x,cluster=factor(clust))
  max.clust <- length(unique(clust))
  mu_k <- matrix(0,ncol=d,nrow=max.clust)
  kappa_k <- numeric(max.clust)
  v_k <- numeric(max.clust)
  T_k <- array(0,c(d,d,max.clust))
  mu <- matrix(0,ncol=d,nrow=max.clust)
  Sigma <- array(0,c(d,d,max.clust))
  for (i in 1:max.clust){
    mu_k[i,] <- (kappa0*mu0+length(clust==i)*colMeans(df[df$cluster==i,][,1:d]))/(kappa0+length(clust==i))
    kappa_k[i] <- kappa0 + length(clust==i)
    v_k[i] <- v0 + length(clust==i)
    y.diff <- sweep(df[df$cluster==i,][,1:d],2,colMeans(df[df$cluster==i,][,1:d]),'-')
    y.diff.sum <- matrix(0,d,d)
    for (j in 1:dim(y.diff)[1]){
      y.diff.sum <- y.diff.sum + t(y.diff[j,])%*%as.matrix(y.diff[j,])
    }
    T_k[,,i] <- T0 + y.diff.sum + kappa0*length(clust==i)/
      (kappa0+length(clust==i))*
      as.matrix((mu0-colMeans(df[df$cluster==i,][,1:d])))%*%t((mu0-colMeans(df[df$cluster==i,][,1:d])))
    Sigma[,,i] <- riwish(v_k[i],T_k[,,i])
    mu[i,] <- rmvnorm(1,mu_k[i,],1/kappa_k[i]*Sigma[,,i])
  }
  n=dim(x)[1]
  ks=max(unique(config_s))
  lpostvalue=numeric(1)
  for (j in 1:ks){
    njs=sum(config_s==j)
    lpostvalue=lpostvalue+sum(dmvnorm(df[df$cluster==j,][,1:d],mu[j,],Sigma[,,j],log=TRUE))+
      log(dwish(Sigma[,,j],v_k[j],solve(T_k[,,j])))+log(dmvnorm(mu[j,],mu_k[j,],1/kappa_k[j]*Sigma[,,j]))
  }   ####dwish will occur NaN, we omit it.
  return(lpostvalue)
}
 




####################sparse mixture model################
####################gamma-normal prior #########
# Hyperparameters   k.true=3   

library(gtools)
K <- 5
d <- 2
g <- 2
a <- 3
delta <- c(rep(0.1,length=K)) # delta < d/2



for (simnum in 1:numofsims){
  begin.time <- proc.time()
  for ( ch in 1:modelchoice){
  if(ch==1){
    niter=2000
    burnin=1000
    result <- gibbs3(x=data3[,,simnum],K=K,d=d,niter=niter,burnin=burnin,a=a,g=g) 
    config=result$z
    post=result$post
  }else{
    burnin=1000
    niter=2000
    dp  <- DirichletProcessMvnormal(scale(data3[,,simnum]))
    result  <- Fit(dp, niter+burnin) 
    config=as.matrix(plyr::ldply(result$labelsChain, rbind))
    post=DPMM.post(data3[,,simnum],result,burnin=burnin,niter=niter)
  }
 



  # Define variables
  
  
  
  # Compute k 
  k=rep(0,niter)
  for(s in 1:niter){
   k[s]=length(unique(config[s,]))
  }
  khatsims[simnum,ch]=mean(k)
  kmodesims[simnum,ch]=c(getmode(k))

######################################################################
#  0-1 loss
########### MAP
  
  

  mapind=which.max(post)
  kMAPsims[simnum,ch]=length(unique(config[mapind,]))
  MAPclus[simnum,,ch]=config[mapind,]
  
  S=dim(config)[1]
  ####remove the empty cluster label###########  
  #remove zero clusters
  for (s in 1:S){
    config_s=config[s,]
    labels_s=unique(config_s)
    for (j in 1:length(labels_s)){
      config[s,config_s==labels_s[j]]=j
    }
  }
  
  config_s=clust.true[simnum,]
  labels_s=unique(config_s)
  for (j in 1:length(labels_s)){
    clust.true[simnum,config_s==labels_s[j]]=j
  }
  

  slpost=sort(post,decreasing=T)
  cb_eps=rep(0,length(alphacb))
  for (ia in 1:length(alphacb)){
    cb_eps[ia]=slpost[ceiling((1-alphacb[ia])*S)]
  }
  cbMAPeps[simnum,,ch]=cb_eps
  MAPsims[simnum,ch]=gibbs.lpost.clusttrue(x=data3[,,simnum],clust=clust.true[simnum,])
  

######################################################################

########### Inference for partitions
psm=comp.psm(config)
output_vi=minVI(psm,config,method=("all"),include.greedy=TRUE,suppress.comment = F,l=n/2)
config_vi=output_vi$cl[1,]
VIclus[simnum,,ch]=config_vi
kVIsims[simnum,ch]=length(unique(config_vi))              ###########k.VI_N^star
VIsims[simnum,ch]=VI(config_vi,matrix(clust.true[simnum,],1,n))     ###########VI(ct,c.star)

cb_vi=credibleball_v2(config_vi, config, c.dist ="VI",alpha=alphacb)
cbVIeps[simnum,,ch]=cb_vi$dist.horiz



output_binder=minbinder.ext(psm,config,method=("all"),include.greedy=TRUE,suppress.comment = F,l=n/2)
config_binder=output_binder$cl[1,]
Binderclus[simnum,,ch]=config_binder
kBindersims[simnum,ch]=length(unique(config_binder))         ###########k.Binder_N^star
Bindersims[simnum,ch]=binder(config_binder,comp.psm(matrix(clust.true[simnum,],1,n)))/(n^2)   ###########Binder(ct,c.star)

cb_binder=credibleball_v2(config_binder, config, c.dist ="Binder",alpha=alphacb)
cbBindereps[simnum,,ch]=cb_binder$dist.horiz

print( paste("Number of simulations completed=", simnum) )
print( paste("k_VI=", kVIsims[simnum,ch]) )
print( paste("k_B=", kBindersims[simnum,ch]) )
print( paste("k_MAP=", kMAPsims[simnum,ch]) )
print( paste("khat=", khatsims[simnum,ch]) )
}
  print(proc.time() - begin.time)
}
#save.image("exampleUNIF_n100.1111.RData")


  

