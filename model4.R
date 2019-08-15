library(dirichletprocess)
library(mvtnorm)
library(gtools)
library(ggplot2)

begin.time <- proc.time()
niter=3000
dp <- DirichletProcessMvnormal(scale(data3))
result6 <- Fit(dp, niter) 
proc.time() - begin.time
plot(result6)


###############plot the weight#################
weightchain <- plyr::ldply(result6$weightsChain, rbind)
weightchain[is.na(weightchain)] <- 0
index <- c(rep(1:dim(weightchain)[1]))
weightchain <- cbind(weightchain,index)
weightchain.melt <- melt(weightchain, id.vars = "index")
colnames(weightchain.melt)[2] <- "components"
ggplot(weightchain.melt , aes(x = index, y = value, color = components)) +
  theme_bw() + ggtitle("Weights in each component") +
  geom_line() +xlab("interation") + ylab("probability")

#################plot the number of component###################

numchain <- plyr::ldply(result6$weightsChain, rbind)
f <- function(x){length(x)-length(x[is.na(x)])}
numchain  <- apply(numchain,1,f)
plot(numchain)

plot(result6) +  ggtitle("Last iteration of the cluster labels")



a <- DPMM.post(x=data3,result=result6)
