
#################mixture 3 univariate normal example###################

library(ggplot2)
library(RColorBrewer)
n=300
aux=runif(n)
clust.true=rep(0,n)
clust.true[aux<0.1]=1
clust.true[aux>0.1& aux<0.7]=2
clust.true[aux>0.7]=3
data4=matrix(0,n,1)
data4[clust.true==1,]=rnorm(sum(clust.true==1),0,1)
data4[clust.true==2,]=rnorm(sum(clust.true==2),5,1)
data4[clust.true==3,]=rnorm(sum(clust.true==3),10,1)

df <- data.frame(mu=data4,cluster=clust.true)

ggplot(df, aes(x=mu)) + geom_density(aes(group=cluster, colour=cluster, fill=cluster), alpha=0.3)





ggplot(df, aes(x=mu, y=c(0))) + 
  geom_point(aes(colour=cluster)) +
  scale_colour_gradientn(colours=brewer.pal(n =3, name = "Pastel2")) + ylim(-0.01,0.5) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),xlim=c(-3,3),colour = "#b3e2cd") +
  stat_function(fun = dnorm, args = list(mean = 5, sd = 1),xlim=c(2,8),colour = "#fdcdac") +
  stat_function(fun = dnorm, args = list(mean = 10, sd = 1),xlim=c(7,13),colour = "#cbd5e8") +
  ylab("density") + ggtitle("Mixture Univariate Normal Data")


data(galaxy)



library(miscF)
y <- as.numeric(data4)
w <- c(0.33, 0.33, 0.34)
mu <- c(rep(mean(data4),length=3))
sigma2 <- c(rep(var(data4),length=3))
Z <- do.call(cbind, lapply(1:3, function(i)
  w[i]*dnorm(y, mu[i], sqrt(sigma2[i]))))
Z <- apply(Z, 1, function(x) which(x==max(x))[1])

result <- uvnm.rjmcmc(y, nsweep=200000, kmax=30, k=3,
                      w, mu, sigma2, Z)

ksave <- result$k.save
round(table(ksave[-(1:100000)])/100000,4)



  focus.k <- 1
  pick.k <- which(ksave==focus.k)
  w <- unlist(result$w.save[pick.k])
  mu <- unlist(result$mu.save[pick.k])
  sigma2 <- unlist(result$sigma2.save[pick.k])
  den.estimate <- rep(w, each=length(y)) *
    dnorm(rep(y, length(w)), mean=rep(mu, each=length(y)),
          sd=rep(sqrt(sigma2), each=length(y)))
  den.estimate <- rowMeans(matrix(den.estimate, nrow=length(y)))*focus.k

  df1 <- data.frame(data4)
  p <- ggplot() + geom_density(df1,aes(x=data4),color="red")
  df <- data.frame(x=sort(y),y=den.estimate[order(y)])
  p <- p+ geom_line(data=df, aes(x=x, y=y), color = "green") 
 
  
  
  focus.k <- 2
  pick.k <- which(ksave==focus.k)
  w <- unlist(result$w.save[pick.k])
  mu <- unlist(result$mu.save[pick.k])
  sigma2 <- unlist(result$sigma2.save[pick.k])
  den.estimate <- rep(w, each=length(y)) *
    dnorm(rep(y, length(w)), mean=rep(mu, each=length(y)),
          sd=rep(sqrt(sigma2), each=length(y)))
  den.estimate <- rowMeans(matrix(den.estimate, nrow=length(y)))*focus.k
  
  df <- data.frame(x=sort(y),y=den.estimate[order(y)])
  p <- p + geom_line(data=df, aes(x=x, y=y), color = "green") 
  
  
   
  for ( i in 3:length(unique(ksave))){
    focus.k <- i
    pick.k <- which(ksave==focus.k)
    w <- unlist(result$w.save[pick.k])
    mu <- unlist(result$mu.save[pick.k])
    sigma2 <- unlist(result$sigma2.save[pick.k])
    den.estimate <- rep(w, each=length(y)) *
      dnorm(rep(y, length(w)), mean=rep(mu, each=length(y)),
            sd=rep(sqrt(sigma2), each=length(y)))
    den.estimate <- rowMeans(matrix(den.estimate, nrow=length(y)))*focus.k
    df <- data.frame(x=sort(y),y=den.estimate[order(y)])
    p <- p + 
    geom_line(data=df,aes(x=x,y=y), color = "midnightblue") 
  
  }

 p + xlab('data')+ylab('denmsity') + ggtitle('Posterior Density with Different Number of Components')
