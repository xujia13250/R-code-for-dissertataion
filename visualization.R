
####################visualization############ Stn_100#################


load("C:/Users/Administrator/Desktop/论文/论文/part3/新建文件夹/new/exampleSTN_n100.final.RData")
library(reshape2)
library(ggplot2)
# Create labels for legend
modelmat=c("model1","model2")
cl=rainbow(modelchoice)

df <- data.frame(sparse=khatsims[,1],DPMM=khatsims[,2])

data<- melt(df)
colnames(data) <- c("model","khat")
ggplot(data,aes(x=khat, fill=model)) + geom_density(alpha=0.25) + ggtitle("Density of khat(Standard normal n=100)")



# Plot distribution MAP value of k
kmin=min(kmodesims)
kmax=max(kmodesims[,1:(modelchoice)])
kheight=matrix(0,modelchoice,kmax-kmin+1)
for (ch in 1:(modelchoice)){
  for(b in 1:(kmax-kmin+1)){
    kheight[ch,b]=sum(kmodesims[,ch]==kmin+b-1)/numofsims
  }
}


kheight <- data.frame(kheight)
kheight$model <- factor(c("Sparse Model","DPMM"), 
                        levels=c("Sparse Model","DPMM"))
colnames(kheight)[1:kmax-kmin+1] <- as.character(rep(1:kmax-kmin+1))
kheight <- melt(kheight, id.vars="model")
ggplot(kheight , aes(variable, value, fill=model)) + 
  geom_bar(stat="identity", position="dodge") + 
  ggtitle("Distribution MAP value of k(Standard normal n=100)") + xlab("estimate k") + ylab('probability')





# Plot distribution VI estimate of k
kmin=min(kVIsims)
kmax=max(max(kVIsims[,1:(modelchoice)]),2)
kheight=matrix(0,modelchoice,kmax-kmin+1)
for (ch in 1:(modelchoice)){
  for(b in 1:(kmax-kmin+1)){
    kheight[ch,b]=sum(kVIsims[,ch]==kmin+b-1)/numofsims
  }
}

kheight <- data.frame(kheight)
kheight$model <- factor(c("Sparse Model","DPMM"), 
                        levels=c("Sparse Model","DPMM"))
colnames(kheight)[1:kmax-kmin+1] <- as.character(rep(1:kmax-kmin+1))
kheight <- melt(kheight, id.vars="model")
ggplot(kheight , aes(variable, value, fill=model)) + 
  geom_bar(stat="identity", position="dodge") + ggtitle("Distribution VI estimate of k(Standard normal n=100)") + xlab("estimate k") + ylab('probability')




# Plot distribution Binder estimate of k
kmin=min(kBindersims)
#kmax=max(kBindersims)
kmax=max(min(20,max(kBindersims)),2)
kheight=matrix(0,modelchoice,kmax-kmin+1)
for (ch in 1:(modelchoice)){
  for(b in 1:(kmax-kmin)){
    kheight[ch,b]=sum(kBindersims[,ch]==kmin+b-1)/numofsims
  }
  kheight[ch,kmax-kmin+1]=sum(kBindersims[,ch]>=kmax)/numofsims
}

kheight <- data.frame(kheight)
kheight$model <- factor(c("Sparse Model","DPMM"), 
                        levels=c("Sparse Model","DPMM"))
colnames(kheight)[1:kmax-kmin+1] <- as.character(rep(1:kmax-kmin+1))
kheight <- melt(kheight, id.vars="model")
ggplot(kheight , aes(variable, value, fill=model)) + 
  geom_bar(stat="identity", position="dodge") + ggtitle("Distribution Binder estimate of k(Standard normal n=100)")+ xlab("estimate k") + ylab('probability')


# Plot distribution MAP estimate of k
kmin=min(kMAPsims)
kmax=max(kMAPsims)
kheight=matrix(0,modelchoice,kmax-kmin+1)
for (ch in 1:(modelchoice)){
  for(b in 1:(kmax-kmin+1)){
    kheight[ch,b]=sum(kMAPsims[,ch]==kmin+b-1)/numofsims
  }
}

kheight <- data.frame(kheight)
kheight$model <- factor(c("Sparse Model","DPMM"), 
                        levels=c("Sparse Model","DPMM"))
colnames(kheight)[1:kmax-kmin+1] <- as.character(rep(1:kmax-kmin+1))
kheight <- melt(kheight, id.vars="model")
ggplot(kheight , aes(variable, value, fill=model)) + 
  geom_bar(stat="identity", position="dodge") + ggtitle("Distribution MAP estimate of k(Standard normal n=100)") + xlab("estimate k") + ylab('probability')



# Compute coverage prob
coverageprob.vi=matrix(0,length(alphacb),modelchoice)
for (a in 1:length(alphacb)){
  for (ch in 1:modelchoice){
    coverageprob.vi[a,ch]=mean(((VIsims[,ch]))<=cbVIeps[,a,ch])
  }
}

coverageprob.binder=matrix(0,length(alphacb),modelchoice)
for (a in 1:length(alphacb)){
  for (ch in 1:modelchoice){
    coverageprob.binder[a,ch]=mean(((Bindersims[,ch]))<=cbBindereps[,a,ch])
  }
}

coverageprob.map=matrix(0,length(alphacb),modelchoice)
for (a in 1:length(alphacb)){
  for (ch in 1:modelchoice){
    coverageprob.map[a,ch]=mean(((MAPsims[,ch]))>=cbMAPeps[,a,ch])
  }
}




df <- data.frame("credible probability"=1-alphacb,"Sparse model"=coverageprob.vi[,1],"DPMM"=coverageprob.vi[,2])
df <- melt(df, id.vars="credible.probability")
colnames(df) <- c("credible.probability","model","coverage.probability")
ggplot(df,aes(x=credible.probability, y=coverage.probability,group=model,color=model)) +
        geom_line()+geom_abline() + ylim(0,1) + ggtitle("Coverage probability of VI(Standard normal n=100)")




df <- data.frame("credible probability"=1-alphacb,"Sparse model"=coverageprob.binder[,1],"DPMM"=coverageprob.binder[,2])
df <- melt(df, id.vars="credible.probability")
colnames(df) <- c("credible.probability","model","coverage.probability")
ggplot(df,aes(x=credible.probability, y=coverage.probability,group=model,color=model)) +
  geom_line()+geom_abline() + ylim(0,1) + ggtitle("Coverage probability of Binder(Standard normal n=100)")



df <- data.frame("credible probability"=1-alphacb,"Sparse model"=coverageprob.map[,1],"DPMM"=coverageprob.map[,2])
df <- melt(df, id.vars="credible.probability")
colnames(df) <- c("credible.probability","model","coverage.probability")
ggplot(df,aes(x=credible.probability, y=coverage.probability,group=model,color=model)) +
  geom_line()+geom_abline() + ylim(0,1) + ggtitle("Coverage probability of MAP(Standard normal n=100)")


df = data.frame(x=data3[,,50],cluster=factor(MAPclus[50,,2]))
ggplot(data=df, aes(x=x.1, y=x.2,color=cluster)) +
  geom_point() + stat_ellipse(level = 0.95, show.legend = F)  + 
  ggtitle("MAP n=100") +xlab("x1") +ylab("x2") 


df = data.frame(x=data3[,,50],cluster=factor(VIclus[50,,2]))
ggplot(data=df, aes(x=x.1, y=x.2,color=cluster)) +
  geom_point() + stat_ellipse(level = 0.95, show.legend = F)  + 
  ggtitle("VI n=100") +xlab("x1") +ylab("x2") 



df = data.frame(x=data3[,,50],cluster=factor(Binderclus[50,,2]))
ggplot(data=df, aes(x=x.1, y=x.2,color=cluster)) +
  geom_point() + stat_ellipse(level = 0.95, show.legend = F)  + 
  ggtitle("Binder n=100") +xlab("x1") +ylab("x2") 





####################visualization############ Stn_300################


load("C:/Users/Administrator/Desktop/论文/论文/part3/新建文件夹/new/exampleSTN_n300.RData")
library(reshape2)
library(ggplot2)
# Create labels for legend
modelmat=c("model1","model2")
cl=rainbow(modelchoice)

df <- data.frame(sparse=khatsims[,1],DPMM=khatsims[,2])

data<- melt(df)
colnames(data) <- c("model","khat")
ggplot(data,aes(x=khat, fill=model)) + geom_density(alpha=0.25) + ggtitle("Density of khat(Standard normal n=300)")



# Plot distribution MAP value of k
kmin=min(kmodesims)
kmax=max(kmodesims[,1:(modelchoice)])
kheight=matrix(0,modelchoice,kmax-kmin+1)
for (ch in 1:(modelchoice)){
  for(b in 1:(kmax-kmin+1)){
    kheight[ch,b]=sum(kmodesims[,ch]==kmin+b-1)/numofsims
  }
}


kheight <- data.frame(kheight)
kheight$model <- factor(c("Sparse Model","DPMM"), 
                        levels=c("Sparse Model","DPMM"))
colnames(kheight)[1:kmax-kmin+1] <- as.character(rep(1:kmax-kmin+1))
kheight <- melt(kheight, id.vars="model")
ggplot(kheight , aes(variable, value, fill=model)) + 
  geom_bar(stat="identity", position="dodge") + ggtitle("Distribution MAP value of k(Standard normal n=300)")+ xlab("estimate k") + ylab('probability')





# Plot distribution VI estimate of k
kmin=min(kVIsims)
kmax=max(max(kVIsims[,1:(modelchoice)]),2)
kheight=matrix(0,modelchoice,kmax-kmin+1)
for (ch in 1:(modelchoice)){
  for(b in 1:(kmax-kmin+1)){
    kheight[ch,b]=sum(kVIsims[,ch]==kmin+b-1)/numofsims
  }
}

kheight <- data.frame(kheight)
kheight$model <- factor(c("Sparse Model","DPMM"), 
                        levels=c("Sparse Model","DPMM"))
colnames(kheight)[1:kmax-kmin+1] <- as.character(rep(1:kmax-kmin+1))
kheight <- melt(kheight, id.vars="model")
ggplot(kheight , aes(variable, value, fill=model)) + 
  geom_bar(stat="identity", position="dodge") + ggtitle("Distribution VI estimate of k(Standard normal n=300)")+ xlab("estimate k") + ylab('probability')




# Plot distribution Binder estimate of k
kmin=min(kBindersims)
#kmax=max(kBindersims)
kmax=max(min(20,max(kBindersims)),2)
kheight=matrix(0,modelchoice,kmax-kmin+1)
for (ch in 1:(modelchoice)){
  for(b in 1:(kmax-kmin)){
    kheight[ch,b]=sum(kBindersims[,ch]==kmin+b-1)/numofsims
  }
  kheight[ch,kmax-kmin+1]=sum(kBindersims[,ch]>=kmax)/numofsims
}

kheight <- data.frame(kheight)
kheight$model <- factor(c("Sparse Model","DPMM"), 
                        levels=c("Sparse Model","DPMM"))
colnames(kheight)[1:kmax-kmin+1] <- as.character(rep(1:kmax-kmin+1))
kheight <- melt(kheight, id.vars="model")
ggplot(kheight , aes(variable, value, fill=model)) + 
  geom_bar(stat="identity", position="dodge") + ggtitle("Distribution Binder estimate of k(Standard normal n=300)")+ xlab("estimate k") + ylab('probability')


# Plot distribution MAP estimate of k
kmin=min(kMAPsims)
kmax=max(kMAPsims)
kheight=matrix(0,modelchoice,kmax-kmin+1)
for (ch in 1:(modelchoice)){
  for(b in 1:(kmax-kmin+1)){
    kheight[ch,b]=sum(kMAPsims[,ch]==kmin+b-1)/numofsims
  }
}

kheight <- data.frame(kheight)
kheight$model <- factor(c("Sparse Model","DPMM"), 
                        levels=c("Sparse Model","DPMM"))
colnames(kheight)[1:kmax-kmin+1] <- as.character(rep(1:kmax-kmin+1))
kheight <- melt(kheight, id.vars="model")
ggplot(kheight , aes(variable, value, fill=model)) + 
  geom_bar(stat="identity", position="dodge") + ggtitle("Distribution MAP estimate of k(Standard normal n=300)")+ xlab("estimate k") + ylab('probability')



# Compute coverage prob
coverageprob.vi=matrix(0,length(alphacb),modelchoice)
for (a in 1:length(alphacb)){
  for (ch in 1:modelchoice){
    coverageprob.vi[a,ch]=mean(((VIsims[,ch]))<=cbVIeps[,a,ch])
  }
}

coverageprob.binder=matrix(0,length(alphacb),modelchoice)
for (a in 1:length(alphacb)){
  for (ch in 1:modelchoice){
    coverageprob.binder[a,ch]=mean(((Bindersims[,ch]))<=cbBindereps[,a,ch])
  }
}

coverageprob.map=matrix(0,length(alphacb),modelchoice)
for (a in 1:length(alphacb)){
  for (ch in 1:modelchoice){
    coverageprob.map[a,ch]=mean(((MAPsims[,ch]))>=cbMAPeps[,a,ch])
  }
}




df <- data.frame("credible probability"=1-alphacb,"Sparse model"=coverageprob.vi[,1],"DPMM"=coverageprob.vi[,2])
df <- melt(df, id.vars="credible.probability")
colnames(df) <- c("credible.probability","model","coverage.probability")
ggplot(df,aes(x=credible.probability, y=coverage.probability,group=model,color=model)) +
  geom_line()+geom_abline() + ylim(0,1) + ggtitle("Coverage probability of VI(Standard normal n=300)")




df <- data.frame("credible probability"=1-alphacb,"Sparse model"=coverageprob.binder[,1],"DPMM"=coverageprob.binder[,2])
df <- melt(df, id.vars="credible.probability")
colnames(df) <- c("credible.probability","model","coverage.probability")
ggplot(df,aes(x=credible.probability, y=coverage.probability,group=model,color=model)) +
  geom_line()+geom_abline() + ylim(0,1) + ggtitle("Coverage probability of Binder(Standard normal n=300)")



df <- data.frame("credible probability"=1-alphacb,"Sparse model"=coverageprob.map[,1],"DPMM"=coverageprob.map[,2])
df <- melt(df, id.vars="credible.probability")
colnames(df) <- c("credible.probability","model","coverage.probability")
ggplot(df,aes(x=credible.probability, y=coverage.probability,group=model,color=model)) +
  geom_line()+geom_abline() + ylim(0,1) + ggtitle("Coverage probability of MAP(Standard normal n=300)")


df = data.frame(x=data3[,,50],cluster=factor(MAPclus[50,,2]))
ggplot(data=df, aes(x=x.1, y=x.2,color=cluster)) +
  geom_point() + stat_ellipse(level = 0.95, show.legend = F)  + 
  ggtitle("MAP n=300") +xlab("x1") +ylab("x2") 


df = data.frame(x=data3[,,50],cluster=factor(VIclus[50,,2]))
ggplot(data=df, aes(x=x.1, y=x.2,color=cluster)) +
  geom_point() + stat_ellipse(level = 0.95, show.legend = F)  + 
  ggtitle("VI n=300") +xlab("x1") +ylab("x2") 



df = data.frame(x=data3[,,50],cluster=factor(Binderclus[50,,2]))
ggplot(data=df, aes(x=x.1, y=x.2,color=cluster)) +
  geom_point() + stat_ellipse(level = 0.95, show.legend = F)  + 
  ggtitle("Binder n=300") +xlab("x1") +ylab("x2") 




####################visualization############ UNIF_100#######################


load("C:/Users/Administrator/Desktop/论文/论文/part3/新建文件夹/new/exampleUNIF_n100.final.RData")
library(reshape2)
library(ggplot2)
# Create labels for legend
modelmat=c("model1","model2")
cl=rainbow(modelchoice)

df <- data.frame(sparse=khatsims[,1],DPMM=khatsims[,2])

data<- melt(df)
colnames(data) <- c("model","khat")
ggplot(data,aes(x=khat, fill=model)) + 
  geom_density(alpha=0.25) + ggtitle("Density of khat(Uniform n=100)")

ggplot(data,aes(x=khat, fill=model)) + 
  geom_density(alpha=0.25) + ggtitle("Density of khat(Uniform n=100)") + xlim(0.5,1.5) 

ggplot(data,aes(x=khat, fill=model)) + 
  geom_density(alpha=0.25) + ggtitle("Density of khat(Uniform n=100)") + xlim(2,6)

# Plot distribution MAP value of k
kmin=min(kmodesims)
kmax=max(kmodesims[,1:(modelchoice)])
kheight=matrix(0,modelchoice,kmax-kmin+1)
for (ch in 1:(modelchoice)){
  for(b in 1:(kmax-kmin+1)){
    kheight[ch,b]=sum(kmodesims[,ch]==kmin+b-1)/numofsims
  }
}


kheight <- data.frame(kheight)
kheight$model <- factor(c("Sparse Model","DPMM"), 
                        levels=c("Sparse Model","DPMM"))
colnames(kheight)[1:kmax-kmin+1] <- as.character(rep(1:kmax-kmin+1))
kheight <- melt(kheight, id.vars="model")
ggplot(kheight , aes(variable, value, fill=model)) + 
  geom_bar(stat="identity", position="dodge") + ggtitle("Distribution MAP value of k(Uniform n=100)")+ xlab("estimate k") + ylab('probability')





# Plot distribution VI estimate of k
kmin=min(kVIsims)
kmax=max(max(kVIsims[,1:(modelchoice)]),2)
kheight=matrix(0,modelchoice,kmax-kmin+1)
for (ch in 1:(modelchoice)){
  for(b in 1:(kmax-kmin+1)){
    kheight[ch,b]=sum(kVIsims[,ch]==kmin+b-1)/numofsims
  }
}

kheight <- data.frame(kheight)
kheight$model <- factor(c("Sparse Model","DPMM"), 
                        levels=c("Sparse Model","DPMM"))
colnames(kheight)[1:kmax-kmin+1] <- as.character(rep(1:kmax-kmin+1))
kheight <- melt(kheight, id.vars="model")
ggplot(kheight , aes(variable, value, fill=model)) + 
  geom_bar(stat="identity", position="dodge") + ggtitle("Distribution VI estimate of k(Uniform n=100)")+ xlab("estimate k") + ylab('probability')




# Plot distribution Binder estimate of k
kmin=min(kBindersims)
#kmax=max(kBindersims)
kmax=max(min(20,max(kBindersims)),2)
kheight=matrix(0,modelchoice,kmax-kmin+1)
for (ch in 1:(modelchoice)){
  for(b in 1:(kmax-kmin)){
    kheight[ch,b]=sum(kBindersims[,ch]==kmin+b-1)/numofsims
  }
  kheight[ch,kmax-kmin+1]=sum(kBindersims[,ch]>=kmax)/numofsims
}

kheight <- data.frame(kheight)
kheight$model <- factor(c("Sparse Model","DPMM"), 
                        levels=c("Sparse Model","DPMM"))
colnames(kheight)[1:kmax-kmin+1] <- as.character(rep(1:kmax-kmin+1))
kheight <- melt(kheight, id.vars="model")
ggplot(kheight , aes(variable, value, fill=model)) + 
  geom_bar(stat="identity", position="dodge") + ggtitle("Distribution Binder estimate of k(Uniform n=100)")+ xlab("estimate k") + ylab('probability')


# Plot distribution MAP estimate of k
kmin=min(kMAPsims)
kmax=max(kMAPsims)
kheight=matrix(0,modelchoice,kmax-kmin+1)
for (ch in 1:(modelchoice)){
  for(b in 1:(kmax-kmin+1)){
    kheight[ch,b]=sum(kMAPsims[,ch]==kmin+b-1)/numofsims
  }
}

kheight <- data.frame(kheight)
kheight$model <- factor(c("Sparse Model","DPMM"), 
                        levels=c("Sparse Model","DPMM"))
colnames(kheight)[1:kmax-kmin+1] <- as.character(rep(1:kmax-kmin+1))
kheight <- melt(kheight, id.vars="model")
ggplot(kheight , aes(variable, value, fill=model)) + 
  geom_bar(stat="identity", position="dodge") + ggtitle("Distribution MAP estimate of k(Uniform n=100)")+ xlab("estimate k") + ylab('probability')



# Compute coverage prob
coverageprob.vi=matrix(0,length(alphacb),modelchoice)
for (a in 1:length(alphacb)){
  for (ch in 1:modelchoice){
    coverageprob.vi[a,ch]=mean(VIsims[,ch]<=cbVIeps[,a,ch])
  }
}

coverageprob.binder=matrix(0,length(alphacb),modelchoice)
for (a in 1:length(alphacb)){
  for (ch in 1:modelchoice){
    coverageprob.binder[a,ch]=mean(((Bindersims[,ch]))<=cbBindereps[,a,ch])
  }
}

coverageprob.map=matrix(0,length(alphacb),modelchoice)
for (a in 1:length(alphacb)){
  for (ch in 1:modelchoice){
    coverageprob.map[a,ch]=mean(((MAPsims[,ch]))>=cbMAPeps[,a,ch])
  }
}

coverageprob.map[is.na(coverageprob.map)]=1


df <- data.frame("credible probability"=1-alphacb,"Sparse model"=coverageprob.vi[,1],"DPMM"=coverageprob.vi[,2])
df <- melt(df, id.vars="credible.probability")
colnames(df) <- c("credible.probability","model","coverage.probability")
ggplot(df,aes(x=credible.probability, y=coverage.probability,group=model,color=model)) +
  geom_line()+geom_abline() + ylim(0,1) + ggtitle("Coverage probability of VI(Uniform n=100)")




df <- data.frame("credible probability"=1-alphacb,"Sparse model"=coverageprob.binder[,1],"DPMM"=coverageprob.binder[,2])
df <- melt(df, id.vars="credible.probability")
colnames(df) <- c("credible.probability","model","coverage.probability")
ggplot(df,aes(x=credible.probability, y=coverage.probability,group=model,color=model)) +
  geom_line()+geom_abline() + ylim(0,1) + ggtitle("Coverage probability of Binder(Uniform n=100)")



df <- data.frame("credible probability"=1-alphacb,"Sparse model"=coverageprob.map[,1],"DPMM"=coverageprob.map[,2])
df <- melt(df, id.vars="credible.probability")
colnames(df) <- c("credible.probability","model","coverage.probability")
ggplot(df,aes(x=credible.probability, y=coverage.probability,group=model,color=model)) +
  geom_line()+geom_abline() + ylim(0,1) + ggtitle("Coverage probability of MAP(Uniform n=100)")



df = data.frame(x=data3[,,50],cluster=factor(MAPclus[50,,2]))
ggplot(data=df, aes(x=x.1, y=x.2,color=cluster)) +
  geom_point() + stat_ellipse(level = 0.95, show.legend = F)  + 
  ggtitle("MAP n=100") +xlab("x1") +ylab("x2") 


df = data.frame(x=data3[,,50],cluster=factor(VIclus[50,,2]))
ggplot(data=df, aes(x=x.1, y=x.2,color=cluster)) +
  geom_point() + stat_ellipse(level = 0.95, show.legend = F)  + 
  ggtitle("VI n=100") +xlab("x1") +ylab("x2") 



df = data.frame(x=data3[,,50],cluster=factor(Binderclus[50,,2]))
ggplot(data=df, aes(x=x.1, y=x.2,color=cluster)) +
  geom_point() + stat_ellipse(level = 0.95, show.legend = F)  + 
  ggtitle("Binder n=100") +xlab("x1") +ylab("x2") 



####################visualization############ UNIF_300################


load("C:/Users/Administrator/Desktop/论文/论文/part3/新建文件夹/new/exampleUNIF_n300.final.RData")
library(reshape2)
library(ggplot2)
# Create labels for legend
modelmat=c("model1","model2")
cl=rainbow(modelchoice)

df <- data.frame(sparse=khatsims[,1],DPMM=khatsims[,2])

data<- melt(df)
colnames(data) <- c("model","khat")
ggplot(data,aes(x=khat, fill=model)) + geom_density(alpha=0.25) + ggtitle("Density of khat(Standard normal n=300)")



# Plot distribution MAP value of k
kmin=min(kmodesims)
kmax=max(kmodesims[,1:(modelchoice)])
kheight=matrix(0,modelchoice,kmax-kmin+1)
for (ch in 1:(modelchoice)){
  for(b in 1:(kmax-kmin+1)){
    kheight[ch,b]=sum(kmodesims[,ch]==kmin+b-1)/numofsims
  }
}


kheight <- data.frame(kheight)
kheight$model <- factor(c("Sparse Model","DPMM"), 
                        levels=c("Sparse Model","DPMM"))
colnames(kheight)[1:kmax-kmin+1] <- as.character(rep(1:kmax-kmin+1))
kheight <- melt(kheight, id.vars="model")
ggplot(kheight , aes(variable, value, fill=model)) + 
  geom_bar(stat="identity", position="dodge") + ggtitle("Distribution MAP value of k(Uniform n=300)")+ xlab("estimate k") + ylab('probability')





# Plot distribution VI estimate of k
kmin=min(kVIsims)
kmax=max(max(kVIsims[,1:(modelchoice)]),2)
kheight=matrix(0,modelchoice,kmax-kmin+1)
for (ch in 1:(modelchoice)){
  for(b in 1:(kmax-kmin+1)){
    kheight[ch,b]=sum(kVIsims[,ch]==kmin+b-1)/numofsims
  }
}

kheight <- data.frame(kheight)
kheight$model <- factor(c("Sparse Model","DPMM"), 
                        levels=c("Sparse Model","DPMM"))
colnames(kheight)[1:kmax-kmin+1] <- as.character(rep(1:kmax-kmin+1))
kheight <- melt(kheight, id.vars="model")
ggplot(kheight , aes(variable, value, fill=model)) + 
  geom_bar(stat="identity", position="dodge") + ggtitle("Distribution VI estimate of k(Uniform n=300)")+ xlab("estimate k") + ylab('probability')




# Plot distribution Binder estimate of k
kmin=min(kBindersims)
#kmax=max(kBindersims)
kmax=max(min(20,max(kBindersims)),2)
kheight=matrix(0,modelchoice,kmax-kmin+1)
for (ch in 1:(modelchoice)){
  for(b in 1:(kmax-kmin)){
    kheight[ch,b]=sum(kBindersims[,ch]==kmin+b-1)/numofsims
  }
  kheight[ch,kmax-kmin+1]=sum(kBindersims[,ch]>=kmax)/numofsims
}

kheight <- data.frame(kheight)
kheight$model <- factor(c("Sparse Model","DPMM"), 
                        levels=c("Sparse Model","DPMM"))
colnames(kheight)[1:kmax-kmin+1] <- as.character(rep(1:kmax-kmin+1))
kheight <- melt(kheight, id.vars="model")
ggplot(kheight , aes(variable, value, fill=model)) + 
  geom_bar(stat="identity", position="dodge") + ggtitle("Distribution Binder estimate of k(Uniform n=300)")+ xlab("estimate k") + ylab('probability')


# Plot distribution MAP estimate of k
kmin=min(kMAPsims)
kmax=max(kMAPsims)
kheight=matrix(0,modelchoice,kmax-kmin+1)
for (ch in 1:(modelchoice)){
  for(b in 1:(kmax-kmin+1)){
    kheight[ch,b]=sum(kMAPsims[,ch]==kmin+b-1)/numofsims
  }
}

kheight <- data.frame(kheight)
kheight$model <- factor(c("Sparse Model","DPMM"), 
                        levels=c("Sparse Model","DPMM"))
colnames(kheight)[1:kmax-kmin+1] <- as.character(rep(1:kmax-kmin+1))
kheight <- melt(kheight, id.vars="model")
ggplot(kheight , aes(variable, value, fill=model)) + 
  geom_bar(stat="identity", position="dodge") + ggtitle("Distribution MAP estimate of k(Uniform n=300)")+ xlab("estimate k") + ylab('probability')



# Compute coverage prob
coverageprob.vi=matrix(0,length(alphacb),modelchoice)
for (a in 1:length(alphacb)){
  for (ch in 1:modelchoice){
    coverageprob.vi[a,ch]=mean(((VIsims[,ch]))<=cbVIeps[,a,ch])
  }
}

coverageprob.binder=matrix(0,length(alphacb),modelchoice)
for (a in 1:length(alphacb)){
  for (ch in 1:modelchoice){
    coverageprob.binder[a,ch]=mean(((Bindersims[,ch]))<=cbBindereps[,a,ch])
  }
}

coverageprob.map=matrix(0,length(alphacb),modelchoice)
for (a in 1:length(alphacb)){
  for (ch in 1:modelchoice){
    coverageprob.map[a,ch]=mean(((MAPsims[,ch]))>=cbMAPeps[,a,ch])
  }
}




df <- data.frame("credible probability"=1-alphacb,"Sparse model"=coverageprob.vi[,1],"DPMM"=coverageprob.vi[,2])
df <- melt(df, id.vars="credible.probability")
colnames(df) <- c("credible.probability","model","coverage.probability")
ggplot(df,aes(x=credible.probability, y=coverage.probability,group=model,color=model)) +
  geom_line()+geom_abline() + ylim(0,1) + ggtitle("Coverage probability of VI(Uniform n=300)")




df <- data.frame("credible probability"=1-alphacb,"Sparse model"=coverageprob.binder[,1],"DPMM"=coverageprob.binder[,2])
df <- melt(df, id.vars="credible.probability")
colnames(df) <- c("credible.probability","model","coverage.probability")
ggplot(df,aes(x=credible.probability, y=coverage.probability,group=model,color=model)) +
  geom_line()+geom_abline() + ylim(0,1) + ggtitle("Coverage probability of Binder(Uniform n=300)")



df <- data.frame("credible probability"=1-alphacb,"Sparse model"=coverageprob.map[,1],"DPMM"=coverageprob.map[,2])
df <- melt(df, id.vars="credible.probability")
colnames(df) <- c("credible.probability","model","coverage.probability")
ggplot(df,aes(x=credible.probability, y=coverage.probability,group=model,color=model)) +
  geom_line()+geom_abline() + ylim(0,1) + ggtitle("Coverage probability of MAP(Uniform n=300)")




df = data.frame(x=data3[,,50],cluster=factor(MAPclus[50,,2]))
ggplot(data=df, aes(x=x.1, y=x.2,color=cluster)) +
  geom_point() + stat_ellipse(level = 0.95, show.legend = F)  + 
  ggtitle("MAP n=300") +xlab("x1") +ylab("x2") 


df = data.frame(x=data3[,,50],cluster=factor(VIclus[50,,2]))
ggplot(data=df, aes(x=x.1, y=x.2,color=cluster)) +
  geom_point() + stat_ellipse(level = 0.95, show.legend = F)  + 
  ggtitle("VI n=300") +xlab("x1") +ylab("x2") 



df = data.frame(x=data3[,,50],cluster=factor(Binderclus[50,,2]))
df[df$cluster==2,]=6
df[df$cluster==5,]=2
df <- df[df$cluster<=5,]
ggplot(data=df, aes(x=x.1, y=x.2,color=cluster)) +
  geom_point() + stat_ellipse(level = 0.95, show.legend = F)  + 
  ggtitle("Binder n=300") +xlab("x1") +ylab("x2") 




####################visualization############ GMM2_100################

load("C:/Users/Administrator/Desktop/论文/论文/part3/新建文件夹/new/exampleGMM2_n100.final.RData")
library(reshape2)
library(ggplot2)
# Create labels for legend
modelmat=c("model1","model2")
cl=rainbow(modelchoice)

df <- data.frame(sparse=khatsims[,1],DPMM=khatsims[,2])

data<- melt(df)
colnames(data) <- c("model","khat")
ggplot(data,aes(x=khat, fill=model)) + geom_density(alpha=0.25) + ggtitle("Density of khat(Guassian Mixture k=2 n=100)")



# Plot distribution MAP value of k
kmin=min(kmodesims)
kmax=max(kmodesims[,1:(modelchoice)])
kheight=matrix(0,modelchoice,kmax-kmin+1)
for (ch in 1:(modelchoice)){
  for(b in 1:(kmax-kmin+1)){
    kheight[ch,b]=sum(kmodesims[,ch]==kmin+b-1)/numofsims
  }
}


kheight <- data.frame(kheight)
kheight$model <- factor(c("Sparse Model","DPMM"), 
                        levels=c("Sparse Model","DPMM"))
colnames(kheight)[1:kmax-kmin+1] <- as.character(rep(1:kmax-kmin+1))
kheight <- melt(kheight, id.vars="model")
ggplot(kheight , aes(variable, value, fill=model)) + 
  geom_bar(stat="identity", position="dodge") + ggtitle("Distribution MAP value of k(Gaussian Mixture k=2 n=100)")+ xlab("estimate k") + ylab('probability')





# Plot distribution VI estimate of k
kmin=min(kVIsims)
kmax=max(max(kVIsims[,1:(modelchoice)]),2)
kheight=matrix(0,modelchoice,kmax-kmin+1)
for (ch in 1:(modelchoice)){
  for(b in 1:(kmax-kmin+1)){
    kheight[ch,b]=sum(kVIsims[,ch]==kmin+b-1)/numofsims
  }
}

kheight <- data.frame(kheight)
kheight$model <- factor(c("Sparse Model","DPMM"), 
                        levels=c("Sparse Model","DPMM"))
colnames(kheight)[1:kmax-kmin+1] <- as.character(rep(1:kmax-kmin+1))
kheight <- melt(kheight, id.vars="model")
ggplot(kheight , aes(variable, value, fill=model)) + 
  geom_bar(stat="identity", position="dodge") + ggtitle("Distribution VI estimate of k(Gaussian Mixture k=2 n=100)")+ xlab("estimate k") + ylab('probability')




# Plot distribution Binder estimate of k
kmin=min(kBindersims)
#kmax=max(kBindersims)
kmax=max(min(20,max(kBindersims)),2)
kheight=matrix(0,modelchoice,kmax-kmin+1)
for (ch in 1:(modelchoice)){
  for(b in 1:(kmax-kmin)){
    kheight[ch,b]=sum(kBindersims[,ch]==kmin+b-1)/numofsims
  }
  kheight[ch,kmax-kmin+1]=sum(kBindersims[,ch]>=kmax)/numofsims
}

kheight <- data.frame(kheight)
kheight$model <- factor(c("Sparse Model","DPMM"), 
                        levels=c("Sparse Model","DPMM"))
colnames(kheight)[1:kmax-kmin+1] <- as.character(rep(1:kmax-kmin+1))
kheight <- melt(kheight, id.vars="model")
ggplot(kheight , aes(variable, value, fill=model)) + 
  geom_bar(stat="identity", position="dodge") + ggtitle("Distribution Binder estimate of k(Gaussian Mixture k=2 n=100)")+ xlab("estimate k") + ylab('probability')


# Plot distribution MAP estimate of k
kmin=min(kMAPsims)
kmax=max(kMAPsims)
kheight=matrix(0,modelchoice,kmax-kmin+1)
for (ch in 1:(modelchoice)){
  for(b in 1:(kmax-kmin+1)){
    kheight[ch,b]=sum(kMAPsims[,ch]==kmin+b-1)/numofsims
  }
}

kheight <- data.frame(kheight)
kheight$model <- factor(c("Sparse Model","DPMM"), 
                        levels=c("Sparse Model","DPMM"))
colnames(kheight)[1:kmax-kmin+1] <- as.character(rep(1:kmax-kmin+1))
kheight <- melt(kheight, id.vars="model")
ggplot(kheight , aes(variable, value, fill=model)) + 
  geom_bar(stat="identity", position="dodge") + ggtitle("Distribution MAP estimate of k(Gaussian Mixture k=2 n=100)")+ xlab("estimate k") + ylab('probability')



# Compute coverage prob
coverageprob.vi=matrix(0,length(alphacb),modelchoice)
for (a in 1:length(alphacb)){
  for (ch in 1:modelchoice){
    coverageprob.vi[a,ch]=mean(((VIsims[,ch]))<=cbVIeps[,a,ch])
  }
}

coverageprob.binder=matrix(0,length(alphacb),modelchoice)
for (a in 1:length(alphacb)){
  for (ch in 1:modelchoice){
    coverageprob.binder[a,ch]=mean(((Bindersims[,ch]))<=cbBindereps[,a,ch])
  }
}

coverageprob.map=matrix(0,length(alphacb),modelchoice)
for (a in 1:length(alphacb)){
  for (ch in 1:modelchoice){
    coverageprob.map[a,ch]=mean(((MAPsims[,ch]))>=cbMAPeps[,a,ch])
  }
}
coverageprob.map[is.na(coverageprob.map)]=1



df <- data.frame("credible probability"=1-alphacb,"Sparse model"=coverageprob.vi[,1],"DPMM"=coverageprob.vi[,2])
df <- melt(df, id.vars="credible.probability")
colnames(df) <- c("credible.probability","model","coverage.probability")
ggplot(df,aes(x=credible.probability, y=coverage.probability,group=model,color=model)) +
  geom_line()+geom_abline() + ylim(0,1) + ggtitle("Coverage probability of VI(Gaussian Mixture k=2 n=100)")




df <- data.frame("credible probability"=1-alphacb,"Sparse model"=coverageprob.binder[,1],"DPMM"=coverageprob.binder[,2])
df <- melt(df, id.vars="credible.probability")
colnames(df) <- c("credible.probability","model","coverage.probability")
ggplot(df,aes(x=credible.probability, y=coverage.probability,group=model,color=model)) +
  geom_line()+geom_abline() + ylim(0,1) + ggtitle("Coverage probability of Binder(Gaussian Mixture k=2 n=100)")



df <- data.frame("credible probability"=1-alphacb,"Sparse model"=coverageprob.map[,1],"DPMM"=coverageprob.map[,2])
df <- melt(df, id.vars="credible.probability")
colnames(df) <- c("credible.probability","model","coverage.probability")
ggplot(df,aes(x=credible.probability, y=coverage.probability,group=model,color=model)) +
  geom_line()+geom_abline() + ylim(0,1) + ggtitle("Coverage probability of MAP(Gaussian Mixture k=2 n=100)")


df = data.frame(x=data3[,,50],cluster=factor(MAPclus[50,,2]))
ggplot(data=df, aes(x=x.1, y=x.2,color=cluster)) +
  geom_point() + stat_ellipse(level = 0.95, show.legend = F)  + 
  ggtitle("MAP n=100") +xlab("x1") +ylab("x2") 


df = data.frame(x=data3[,,50],cluster=factor(VIclus[50,,2]))
ggplot(data=df, aes(x=x.1, y=x.2,color=cluster)) +
  geom_point() + stat_ellipse(level = 0.95, show.legend = F)  + 
  ggtitle("VI n=100") +xlab("x1") +ylab("x2") 



df = data.frame(x=data3[,,50],cluster=factor(Binderclus[50,,2]))
ggplot(data=df, aes(x=x.1, y=x.2,color=cluster)) +
  geom_point() + stat_ellipse(level = 0.95, show.legend = F)  + 
  ggtitle("Binder n=100") +xlab("x1") +ylab("x2") 



df = data.frame(x=data3[,,50],cluster=factor(clust.true[50,]))
ggplot(data=df, aes(x=x.1, y=x.2,color=cluster)) +
  geom_point() + stat_ellipse(level = 0.95, show.legend = F)  + 
  ggtitle("MAP n=100") +xlab("x1") +ylab("x2") 




####################visualization############ GMM2_300###############

load("C:/Users/Administrator/Desktop/论文/论文/part3/新建文件夹/new/exampleGMM2_n300.RData")
library(reshape2)
library(ggplot2)
# Create labels for legend
modelmat=c("model1","model2")
cl=rainbow(modelchoice)

df <- data.frame(sparse=khatsims[,1],DPMM=khatsims[,2])

data<- melt(df)
colnames(data) <- c("model","khat")
ggplot(data,aes(x=khat, fill=model)) + geom_density(alpha=0.25) + ggtitle("Density of khat(Guassian Mixture k=2 n=300)")



# Plot distribution MAP value of k
kmin=min(kmodesims)
kmax=max(kmodesims[,1:(modelchoice)])
kheight=matrix(0,modelchoice,kmax-kmin+1)
for (ch in 1:(modelchoice)){
  for(b in 1:(kmax-kmin+1)){
    kheight[ch,b]=sum(kmodesims[,ch]==kmin+b-1)/numofsims
  }
}


kheight <- data.frame(kheight)
kheight$model <- factor(c("Sparse Model","DPMM"), 
                        levels=c("Sparse Model","DPMM"))
colnames(kheight)[1:kmax-kmin+1] <- as.character(rep(1:kmax-kmin+1))
kheight <- melt(kheight, id.vars="model")
ggplot(kheight , aes(variable, value, fill=model)) + 
  geom_bar(stat="identity", position="dodge") + ggtitle("Distribution MAP value of k(Gaussian Mixture k=2 n=300)")+ xlab("estimate k") + ylab('probability')





# Plot distribution VI estimate of k
kmin=min(kVIsims)
kmax=max(max(kVIsims[,1:(modelchoice)]),2)
kheight=matrix(0,modelchoice,kmax-kmin+1)
for (ch in 1:(modelchoice)){
  for(b in 1:(kmax-kmin+1)){
    kheight[ch,b]=sum(kVIsims[,ch]==kmin+b-1)/numofsims
  }
}

kheight <- data.frame(kheight)
kheight$model <- factor(c("Sparse Model","DPMM"), 
                        levels=c("Sparse Model","DPMM"))
colnames(kheight)[1:kmax-kmin+1] <- as.character(rep(1:kmax-kmin+1))
kheight <- melt(kheight, id.vars="model")
ggplot(kheight , aes(variable, value, fill=model)) + 
  geom_bar(stat="identity", position="dodge") + ggtitle("Distribution VI estimate of k(Gaussian Mixture k=2 n=300)")+ xlab("estimate k") + ylab('probability')




# Plot distribution Binder estimate of k
kmin=min(kBindersims)
#kmax=max(kBindersims)
kmax=max(min(20,max(kBindersims)),2)
kheight=matrix(0,modelchoice,kmax-kmin+1)
for (ch in 1:(modelchoice)){
  for(b in 1:(kmax-kmin)){
    kheight[ch,b]=sum(kBindersims[,ch]==kmin+b-1)/numofsims
  }
  kheight[ch,kmax-kmin+1]=sum(kBindersims[,ch]>=kmax)/numofsims
}

kheight <- data.frame(kheight)
kheight$model <- factor(c("Sparse Model","DPMM"), 
                        levels=c("Sparse Model","DPMM"))
colnames(kheight)[1:kmax-kmin+1] <- as.character(rep(1:kmax-kmin+1))
kheight <- melt(kheight, id.vars="model")
ggplot(kheight , aes(variable, value, fill=model)) + 
  geom_bar(stat="identity", position="dodge") + ggtitle("Distribution Binder estimate of k(Gaussian Mixture k=2 n=300)")+ xlab("estimate k") + ylab('probability')


# Plot distribution MAP estimate of k
kmin=min(kMAPsims)
kmax=max(kMAPsims)
kheight=matrix(0,modelchoice,kmax-kmin+1)
for (ch in 1:(modelchoice)){
  for(b in 1:(kmax-kmin+1)){
    kheight[ch,b]=sum(kMAPsims[,ch]==kmin+b-1)/numofsims
  }
}

kheight <- data.frame(kheight)
kheight$model <- factor(c("Sparse Model","DPMM"), 
                        levels=c("Sparse Model","DPMM"))
colnames(kheight)[1:kmax-kmin+1] <- as.character(rep(1:kmax-kmin+1))
kheight <- melt(kheight, id.vars="model")
ggplot(kheight , aes(variable, value, fill=model)) + 
  geom_bar(stat="identity", position="dodge") + ggtitle("Distribution MAP estimate of k(Gaussian Mixture k=2 n=300)")+ xlab("estimate k") + ylab('probability')



# Compute coverage prob
coverageprob.vi=matrix(0,length(alphacb),modelchoice)
for (a in 1:length(alphacb)){
  for (ch in 1:modelchoice){
    coverageprob.vi[a,ch]=mean(((VIsims[,ch]))<=cbVIeps[,a,ch])
  }
}

coverageprob.binder=matrix(0,length(alphacb),modelchoice)
for (a in 1:length(alphacb)){
  for (ch in 1:modelchoice){
    coverageprob.binder[a,ch]=mean(((Bindersims[,ch]))<=cbBindereps[,a,ch])
  }
}

coverageprob.map=matrix(0,length(alphacb),modelchoice)
for (a in 1:length(alphacb)){
  for (ch in 1:modelchoice){
    coverageprob.map[a,ch]=mean(((MAPsims[,ch]))>=cbMAPeps[,a,ch])
  }
}




df <- data.frame("credible probability"=1-alphacb,"Sparse model"=coverageprob.vi[,1],"DPMM"=coverageprob.vi[,2])
df <- melt(df, id.vars="credible.probability")
colnames(df) <- c("credible.probability","model","coverage.probability")
ggplot(df,aes(x=credible.probability, y=coverage.probability,group=model,color=model)) +
  geom_line()+geom_abline() + ylim(0,1) + ggtitle("Coverage probability of VI(Gaussian Mixture k=2 n=300)")




df <- data.frame("credible probability"=1-alphacb,"Sparse model"=coverageprob.binder[,1],"DPMM"=coverageprob.binder[,2])
df <- melt(df, id.vars="credible.probability")
colnames(df) <- c("credible.probability","model","coverage.probability")
ggplot(df,aes(x=credible.probability, y=coverage.probability,group=model,color=model)) +
  geom_line()+geom_abline() + ylim(0,1) + ggtitle("Coverage probability of Binder(Gaussian Mixture k=2 n=300)")



df <- data.frame("credible probability"=1-alphacb,"Sparse model"=coverageprob.map[,1],"DPMM"=coverageprob.map[,2])
df <- melt(df, id.vars="credible.probability")
colnames(df) <- c("credible.probability","model","coverage.probability")
ggplot(df,aes(x=credible.probability, y=coverage.probability,group=model,color=model)) +
  geom_line()+geom_abline() + ylim(0,1) + ggtitle("Coverage probability of MAP(Gaussian Mixture k=2 n=300)")


df = data.frame(x=data3[,,50],cluster=factor(MAPclus[50,,2]))
ggplot(data=df, aes(x=x.1, y=x.2,color=cluster)) +
  geom_point() + stat_ellipse(level = 0.95, show.legend = F)  + 
  ggtitle("MAP n=300") +xlab("x1") +ylab("x2") 


df = data.frame(x=data3[,,50],cluster=factor(VIclus[50,,2]))
ggplot(data=df, aes(x=x.1, y=x.2,color=cluster)) +
  geom_point() + stat_ellipse(level = 0.95, show.legend = F)  + 
  ggtitle("VI n=300") +xlab("x1") +ylab("x2") 



df = data.frame(x=data3[,,50],cluster=factor(Binderclus[50,,2]))
ggplot(data=df, aes(x=x.1, y=x.2,color=cluster)) +
  geom_point() + stat_ellipse(level = 0.95, show.legend = F)  + 
  ggtitle("Binder n=300") +xlab("x1") +ylab("x2") 




