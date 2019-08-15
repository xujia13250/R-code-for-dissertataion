gibbs3 = function(x,K,d,niter,a,g,burnin){
  
  ##function ###
  #############################################################
  
  
  
  normalize = function(x){return(x/sum(x))}
  
  #' @param x an n vector of data
  #' @param eta a k vector
  #' @param mu a k vector
  #' @param Sigmamns1 a d*(K*d) matrix
  
  sample_z = function(x,mu,Sigmamns1,eta){
    z.pro <- matrix(NA,ncol=K,nrow=dim(x)[1])
    for (i in 1:K){
      z.pro[,i]<-dmvnorm(x,mu[i,],solve(Sigmamns1[,((i-1)*d+1):(i*d)]),log=TRUE)
    }
    p.zz.given.x <- sweep(z.pro,2,log(as.vector(eta)),'+')
    p.zz.given.x.max <- apply(p.zz.given.x, 1, function(x) max(x, na.rm = TRUE))
    p.z.given.x <- exp(sweep(p.zz.given.x,1,p.zz.given.x.max,'-'))
    p.z.given.x <- t(apply(p.z.given.x,1,normalize))
    z = rep(0, dim(x)[1])
    for(i in 1:length(z)){
      z[i] = sample(1:length(eta), size=1,prob=p.z.given.x[i,],replace=TRUE)
    }
    return(z)
  }
  
  #############################################################
  #' @param g,a,h scalar 
  #' @param K the number o clusters
  #' @param d the dimension
  #' @param mu the K*d matrix 
  #' @param Sigmamns1 the d*(K*d) matrix
  sample_G0 = function(g, K, d, a, h,Sigmamns1){
    Sigmamns1.sum <-matrix(0,nrow=d,ncol=d)
    dim(Sigmamns1)=c(d,d,K)
    for ( i in 1:K){
      Sigmamns1.sum <- Sigmamns1.sum + Sigmamns1[,,i]
    }
    G0 <- rWishart(1,2*g+2*K*a,solve(2*h+2*Sigmamns1.sum))[,,1]
    return(G0)
  }
  
  #############################################################
  
  #' @param z an n vector of cluster allocations (1...k)
  #' @param K the number of clusters
  #' @param delta the parameters in Dirichelet
  sample_eta = function(z,K,delta){
    counts = colSums(outer(z,1:K,FUN="=="))
    eta = gtools::rdirichlet(1,counts+delta)
    return(eta)
  }
  #############################################################
  
  
  #' @param x an n vector of data
  #' @param z an n vector of cluster allocations
  #' @param K the number o clusters
  #' @param omega the prior mean for mu
  #' @param Omega the prior var for mu
  #' @param d the dimension
  #' @param Sigmamns1 the d*(K*d) matrix 
  sample_mu = function(x, z, K,d, omega,Omega, Sigmamns1){
    mu <- matrix(NA,nrow=K,ncol=d)
    df = data.frame(x=x,z=z)
    for(i in 1:K){
      sample.size = sum(z==i)
      if(sample.size==0){
        sample.sum = matrix(0,nrow=1,ncol=d)
      }else{
        sample.sum=t(as.matrix(apply(df[df$z==i,],2,sum)[1:d]))
      }
      post.var = solve(sample.size*Sigmamns1[,((i-1)*d+1):(i*d)]+Omega)
      post.mean = (omega %*% t(Omega) + sample.sum%*%t(Sigmamns1[,((i-1)*d+1):(i*d)]))%*%t(solve(sample.size*Sigmamns1[,((i-1)*d+1):(i*d)]+Omega))
      mu[i,] <- rmvnorm(1,post.mean, post.var)
    }
    return(mu)
  }
  #############################################################
  
  #' @param x an n vector of data
  #' @param z an n vector of cluster allocations
  #' @param K the number o clusters
  #' @param d the dimension
  #' @param mu the K*d matrix 
  #' @param G0 the d*d matrix of Sigma prior
  sample_Sigmamns1 = function(x, z, K,d,G0,mu){
    df = data.frame(x=x,z=z)
    sample.size = sum(z==1)
    if(sample.size==0){
      sample.diff.sum.mat = matrix(0,nrow=d,ncol=d)
    }else{
      sample.diff.sum.mat = matrix(0,nrow=d,ncol=d)
      sample.diff <- sweep(df[z==1,][,1:d],2,mu[1,],'-')
      for ( j in 1:dim(sample.diff)[1]){
        sample.diff.sum.mat = sample.diff.sum.mat + t(sample.diff[j,])%*%as.matrix(sample.diff[j,])
      }
    }
    Sigmamns1 <- rWishart(1,2*a+sample.size,solve(2*G0+ sample.diff.sum.mat))[,,1]
    for(i in 2:K){
      sample.size = sum(z==i)
      if(sample.size==0){
        sample.diff.sum.mat = matrix(0,nrow=d,ncol=d)
      }else{
        sample.diff.sum.mat = matrix(0,nrow=d,ncol=d)
        sample.diff <- sweep(df[z==i,][,1:d],2,mu[i,],'-')
        for ( j in 1:dim(sample.diff)[1]){
          sample.diff.sum.mat = sample.diff.sum.mat + t(sample.diff[j,])%*%as.matrix(sample.diff[j,])
        }
      }
      Sigmamns1 <- cbind(Sigmamns1,rWishart(1,2*a+sample.size,solve(2*G0+ sample.diff.sum.mat))[,,1])
    }
    return(Sigmamns1)
  }
  
  ###################loglik#######################
  loglik <- function(x,mu,Sigmamns1,z){
    d <- dim(x)[2]
    tem.log <- numeric(1)
    df <- data.frame(x=x,z=z)
    for ( i in 1:max(df$z)){
      tem.log <- tem.log + sum(dmvnorm(df[df$z==i,][,1:d],mean=mu[i,],sigma=solve(Sigmamns1[,((i-1)*d+1):(i*d)]),log = TRUE))
    }
    return(tem.log)
  }
  ######################AIC#####################
  
  AIC <- function(x,mu,Sigmamns1,z){
    d <- dim(x)[2]
    tem.log <- numeric(1)
    df <- data.frame(x=x,z=z)
    for ( i in 1:max(df$z)){
      tem.log <- tem.log + sum(dmvnorm(df[df$z==i,][,1:d],mean=mu[i,],sigma=solve(Sigmamns1[,((i-1)*d+1):(i*d)]),log = TRUE))
    }
    return(-2*tem.log+2*2*3)
  }
  
  ######################BIC################################
  BIC <- function(x,mu,Sigmamns1,z){
    d <- dim(x)[2]
    n <- dim(x)[1]
    tem.log <- numeric(1)
    df <- data.frame(x=x,z=z)
    for ( i in 1:max(df$z)){
      tem.log <- tem.log + sum(dmvnorm(df[df$z==i,][,1:d],mean=mu[i,],sigma=solve(Sigmamns1[,((i-1)*d+1):(i*d)]),log = TRUE))
    }
    return(-2*tem.log+2*3*log(n))
  }
  
  #####################post###################################
  post <- function(loglik,eta){
    post <- loglik + sum(log(eta))
  }
  
  
  
  ##### run the gibbs##############
  # initial setting
  ##allocate the data to each component
  eta <- rdirichlet(1, delta)
  n <- dim(x)[1]
  z <- apply(matrix(runif(n),nrow=n,ncol=(K-1))>t(matrix(eta[1:(K-1)],nrow=(K-1),ncol=n)),1,sum)+1
  
  ## mu and var
  
  m0 <- apply(x,2,median)
  M0 <- diag(100,d,d)
  lambda.d <- rgamma(d,1,1)
  Gamma.p <- diag(sqrt(lambda.d),d,d)
  R0 <- diag(c((max(x[,1]) - min(x[,1]))^(-2),(max(x[,2]) - min(x[,2]))^(-2)),nrow=d)
  Omega <- Gamma.p%*%R0%*%Gamma.p
  omega <- rmvnorm(1,m0,M0)
  mu <- rmvnorm(K,omega,Omega)
  h <- 100*g*R0/a
  G0 <- rWishart(1,2*g,solve(2*h))[,,1]
  Sigmamns1<- rWishart(1,2*a,solve(2*G0))[,,1]

  for (i in 1:(K-1)){
    Sigmamns1 <- cbind(Sigmamns1,rWishart(1,2*a,solve(2*G0))[,,1])
  }
  
  AIC.DATA <- AIC(x=x,mu=mu,Sigmamns1 = Sigmamns1, z=z)
  BIC.DATA <- BIC(x=x,mu=mu,Sigmamns1 = Sigmamns1, z=z)
  loglikhood.data <- loglik(x=x,mu=mu,Sigmamns1 = Sigmamns1, z=z)
  post.data <- post(loglik=loglikhood.data,eta=eta)
  
  
  res = list(mu=matrix(NA,nrow=niter+burnin, ncol=K*d), eta = matrix(NA,nrow=niter+burnin,ncol=K), 
             z = matrix(NA,nrow=niter+burnin, ncol=dim(x)[1]),Sigmamns1=matrix(NA,nrow=(niter+burnin)*d, ncol=K*d),
             G0=matrix(NA,nrow=(niter+burnin)*d, ncol=d), loglikhood=numeric(niter+burnin),
             AIC=numeric(niter+burnin),BIC=numeric(niter+burnin),post=numeric(niter+burnin))
  
  
  
  for ( i in 1:K){
    res$mu[1,((i-1)*d+1):(i*d)]=mu[i,]
  }
  res$eta[1,]=eta
  res$z[1,]=z 
  res$Sigmamns1[1:d,]=Sigmamns1
  res$G0[1:d,]=G0
  res$AIC[1]=AIC.DATA
  res$BIC[1]=BIC.DATA
  res$loglikhood[1]=loglikhood.data
  res$post[1]=post.data
  
  
  for(i in 2:(burnin+niter)){
    z = sample_z (x=x,mu=mu,Sigmamns1=Sigmamns1,eta=eta)
    G0=sample_G0(g=g, K=K,d=d, a=a,h=h,Sigmamns1=Sigmamns1)
    eta = sample_eta(z=z,K=K,delta=delta)
    mu = sample_mu(x=x, z=z, K=K,d=d, omega=omega,Omega=Omega, Sigmamns1=Sigmamns1)
    Sigmamns1=sample_Sigmamns1(x=x, z=z, K=K,d=d,mu=mu,G0=G0)
    AIC.DATA <- AIC(x=x,mu=mu,Sigmamns1 = Sigmamns1, z=z)
    BIC.DATA <- BIC(x=x,mu=mu,Sigmamns1 = Sigmamns1, z=z)
    loglikhood.data <- loglik(x=x,mu=mu,Sigmamns1 = Sigmamns1, z=z)
    post.data <- post(loglik=loglikhood.data,eta=eta)
    for ( j in 1:K){
      res$mu[i,((j-1)*d+1):(j*d)]=mu[j,]
    }
    res$eta[i,] = eta
    res$z[i,] = z
    res$Sigmamns1[((i-1)*d+1):(i*d),] = Sigmamns1
    res$G0[((i-1)*d+1):(i*d),] = G0
    res$AIC[i]=AIC.DATA
    res$BIC[i]=BIC.DATA
    res$loglikhood[i]=loglikhood.data 
    res$post[i]=post.data
  }
  res$mu <- res$mu[-c(1:burnin),]
  res$eta <- res$eta[-c(1:burnin),]
  res$z <- res$z[-c(1:burnin),]
  res$G0 <- res$G0[-c(1:burnin),]
  res$Sigmamns1 <- res$Sigmamns1[-c(1:burnin*d),]
  res$AIC<- res$AIC[-c(1:burnin)]
  res$BIC <- res$BIC[-c(1:burnin)]
  res$loglikhood<-res$loglikhood[-c(1:burnin)]
  res$post <- res$post[-c(1:burnin)]
  return(res)
}


####################gamma-normal prior #########
# Hyperparameters   k.true=3   
library(gtools)
K <- 7
d <- 2
g <- 2
a <- 3
delta <- c(rep(0.1,length=K)) # delta < d/2

niter=3000
burnin=2000
begin.time <- proc.time()
result4 <- gibbs3(x=data3,K=K,d=d,niter=niter,burnin=burnin,a=a,g=g) 
proc.time() - begin.time




###############standard normal prior ##########
library(gtools)
K <- 7
d <- 2
g <- 2
a <- 3
delta <- c(rep(0.1,length=K)) # delta < d/2

niter=5000
burnin=0
begin.time <- proc.time()
result5 <- gibbs1(x=data3,K=K,d=d,niter=niter,burnin=burnin,a=a,g=g) 
proc.time() - begin.time

################plot estimate mu#########################
library(ggplot2)
library(reshape2)
library(RColorBrewer)
mu.trans <- function(x,d){
  n <- dim(x)[1]
  K <- dim(x)[2]/d
  d1 <- x[,1:d]
  for (i in 2:K){
    d1 <- rbind(d1,x[,((i-1)*d+1):(i*d)])
  }
  d2 <- data.frame(x=d1, z=rep(1:K,each=n))
  return(d2)
}
#####################gamma-normal prior###################

#plot the mu , delete the burn-in period 
mu.df.1 <- mu.trans(result4$mu[-c(1:1999),],2)##seperate the true cluster
mu.df.1[which(mu.df.1[,3]==1),][,3] <-8
mu.df.1[which(mu.df.1[,3]==2),][,3] <-1
mu.df.1[which(mu.df.1[,3]==4),][,3] <-1
mu.df.1[which(mu.df.1[,3]==7),][,3] <-1
mu.df.1[which(mu.df.1[,3]==6),][,3] <-1
mu.df.1[which(mu.df.1[,3]==8),][,3] <-2
mu.df.1[which(mu.df.1[,3]==5),][,3] <-4
mu.df.1 <- mu.df.1[order(mu.df.1[,3]),]  #change the order
ggplot(mu.df.1,aes(x=mu.df.1$x.1,y=mu.df.1$x.2)) + geom_point(aes(colour=z,alpha=z)) +
  scale_colour_gradientn(colours=brewer.pal(n = 4, name = "Set3")) + 
  theme_dark() +
  ggtitle("Estimate mu with gamma-normal prior") + xlab("mu1") +ylab("mu2") 

df4 <- data.frame(index=rep(1:(niter-1)),
                  prob=result4$eta)
didfm3 <- melt(df4,id.vars="index")
ggplot(didfm3,aes(x=index,y=value)) + geom_line(aes(color=variable)) + ggtitle("trace plot for weights with gamma-normal prior")

#####################standard normal prior################

#plot the mu , delete the burn-in period   
mu.df.1 <- mu.trans(result5$mu[-c(1:1999),],2)
mu.df.1[which(mu.df.1[,3]==1),][,3]  <-8
mu.df.1[which(mu.df.1[,3]==2),][,3]  <-1
mu.df.1[which(mu.df.1[,3]==5),][,3]  <-1
mu.df.1[which(mu.df.1[,3]==6),][,3]  <-1
mu.df.1[which(mu.df.1[,3]==7),][,3]  <-1
mu.df.1[which(mu.df.1[,3]==8),][,3]  <-2
mu.df.1 <- mu.df.1[order(mu.df.1[,3]),]  #change the order
ggplot(mu.df.1,aes(x=mu.df.1$x.1,y=mu.df.1$x.2)) + geom_point(aes(colour=z,alpha=z)) +
  scale_colour_gradientn(colours=brewer.pal(n = 4, name = "Set3")) + 
  ggtitle("Estimate mu with standard normal prior") + xlab("mu1") +ylab("mu2") +
  theme_dark() 


df5 <- data.frame(index=rep(1:(niter-1)),
                  prob=result5$eta)
didfm5 <- melt(df5,id.vars="index")
ggplot(didfm5,aes(x=index,y=value)) + geom_line(aes(color=variable)) + ggtitle("trace plot for weights with standard normal prior")

############plot each chain of components###############
####standard normal prior##########
library(reshape2)
library(ggplot2)
a <- apply(result5$z,1,table)
c.chain <- plyr::ldply(a, rbind)
c.chain[is.na(c.chain)] <- 0
index <- c(rep(1:dim(c.chain)[1]))
c.chain <- cbind(c.chain,index)
c.chain.melt <- melt(c.chain, id.vars = "index")
colnames(c.chain.melt)[2] <- "components"
ggplot(c.chain.melt , aes(x = index, y = value, color = components)) +
  theme_bw() + ggtitle("Total number in each component(standard normal prior)") +
  geom_line() +xlab("interation") + ylab("number in components") + labs(fill ="components")

#####gamma-normal prior##############
library(reshape2)
library(ggplot2)
a <- apply(result4$z,1,table)
c.chain <- plyr::ldply(a, rbind)
c.chain[is.na(c.chain)] <- 0
index <- c(rep(1:dim(c.chain)[1]))
c.chain <- cbind(c.chain,index)
c.chain.melt <- melt(c.chain, id.vars = "index")
colnames(c.chain.melt)[2] <- "components"
ggplot(c.chain.melt , aes(x = index, y = value, color = components)) +
  theme_bw() + ggtitle("Total number in each component(gamma-normal prior)") +
  geom_line() +xlab("interation") + ylab("number in components") 

