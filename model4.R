### Function that obtains posterior samples from  DP mixture model

###INPUT
# S:number of iterations to save
# burnin: number of iterations to discard
# x: observed covariates (nxp)
# mu_0, c_x: prior parameters for mu  (px1 and px1)
# a_x, b_x: prior parameters for sigma^2_x (px1 and px1)
# u_alpha, v_alpha: prior parameters for alpha (1x1 and 1x1)
# config_init: initial configuration (1xn)
# mu_x_init, sigma_x_init: initial values of mu_x and sigma_x for each y-x cluster (list with kn_y elements, each element is pxkn_x[j])
# alpha_init: initial value of mass parameter (1x1)

###OUPUT
# config: sampled configurations (Sxn)
# mu_x, sigma_x: sampled x cluster parameters (list with S elements, each element is pxkn[s] and pxkn[s])
# alpha: sampled precision parameters (Sx1)


dp_mcmc=function(S, burnin, x, mu_0, c_x, a_x, b_x, alpha, config_init, mu_x_init, sigma_x_init){

#define n and p
n=dim(x)[1]
p=dim(x)[2]

#Empty matrices to store MCMC output

config=matrix(0,S,n)
mu_x=list()
sigma_x=list()

#Initialize
config_s=config_init
mu_x_s=mu_x_init
sigma_x_s=sigma_x_init

#terms that don't depend on s
mucmu_x=c_x*mu_0^2

for(s in 1:(S+burnin-1)){

	#Sample configurations at s+1

	#configuration to be updated
	c_update=config_s

	#unique parameters to be updated if new value is chosen
	sigmax_update=sigma_x_s
	mux_update=mu_x_s

	#updated hyperparameters that dont depend on i	
	a_x_hat=a_x+1/2
	c_x_hat=c_x+1

	for(i in 1:n){
	
		#Define configuration without i
		c_ni=c_update[-i]
		c_i=c_update[i]
		uniq_c_ni=unique(c_ni)
		dn_ni=length(uniq_c_ni)
	
		#check if c_i is a class with the single element i
		singleton=(sum(c_ni==c_i)==0)


		if(singleton){
			#remove the parameters for class c_i
			sigmax_update[,c_i]=sigmax_update[,dn_ni+1]
			sigmax_update=matrix(sigmax_update[,-(dn_ni+1)], nrow=p)
			mux_update[,c_i]=mux_update[,dn_ni+1]
			mux_update=matrix(mux_update[,-(dn_ni+1)],nrow=p)
			#change c_ni, updated config
			c_ni[c_ni==(dn_ni+1)]=c_i
			c_update[c_update==(dn_ni+1)]=c_i
			c_update[i]=dn_ni+1		
		}


		#Set value of config for observation i
		log_prob=matrix(0,dn_ni+1,1)
		log_x=matrix(0,dn_ni+1,1)
		mux_update=matrix(mux_update, nrow=p, ncol=dn_ni)
		sigmax_update=matrix(sigmax_update, nrow=p, ncol=dn_ni)


		#calculate probabilities for old values
		for(j in 1:dn_ni){
			n_j_ni=sum(c_ni==j)
			log_x[j]=sum(dnorm(x[i,], mux_update[,j], sigmax_update[,j]^.5, log=TRUE))
			log_prob[j]=log(n_j_ni/(alpha+n-1))+log_x[j]
		}

		#calculate probability of a new value
		#updated parameters
		x_sc=(x[i,]-mu_0)*(c_x*a_x/(c_x_hat*b_x))^.5
		log_x[dn_ni+1]=sum(dt(x_sc,2*a_x, log=TRUE)+.5*(log(c_x)-log(c_x+1)-log(b_x)+log(a_x)))
		log_prob[dn_ni+1]=log(alpha/(alpha+n-1))+log_x[dn_ni+1]
	
		#need to find normalzing constant
		#f <- function (k) sum(exp(k+log_prob))-1
		#k_star=uniroot(f,lower=-10^(10),upper=10^20)$root
		#probs=exp(k_star+log_prob)
		probs=exp(log_prob)/sum(exp(log_prob))

		#set value of config according to probability
 		c_update[i]=sum(runif(1)>c(cumsum(probs[-(dn_ni+1)]),1))+1

		
		#update parameters if new chosen
		if(c_update[i]==(dn_ni+1)){
			#updated parameters
			mu_x_hat=1/c_x_hat*(c_x*mu_0+x[i,])
			b_x_hat=b_x+(x[i,]^2+mucmu_x-c_x_hat*mu_x_hat^2)/2

			sigmax_new=rinvgamma(p,a_x_hat, b_x_hat)
			mu_x_new=rnorm(p, mu_x_hat, (sigmax_new/c_x_hat)^.5)
			mux_update=cbind(mux_update,mu_x_new)
			sigmax_update=cbind(sigmax_update,sigmax_new)
		}	
	}


	config_s=c_update
	dn=length(unique(c_update))

	#Draw new parameters for each cluster
	mu_x_s=matrix(0,p,dn)
	sigma_x_s=matrix(0,p,dn)
	for(c in 1:dn){

		#observations in class c
		n_c=sum(config_s==c)
		x_c=matrix(x[config_s==c],nrow=n_c)
	
		c_x_hat= c_x+n_c
		mux_hat= 1/c_x_hat*(c_x*mu_0+ colSums(x_c))
		
		#for sigma_x
		a_x_hat=a_x+n_c/2
		b_x_hat=b_x+(colSums(x_c^2)+mucmu_x-c_x_hat*mux_hat^2)/2

		#draw value from posterior
		sigma_x_s[,c]=rinvgamma(p,a_x_hat, b_x_hat)
		mu_x_s[,c]=rnorm(p,mux_hat, (sigma_x_s[,c]/c_x_hat)^.5)

	}

	###Draw other parameters

	### alpha


	if (s%%1000==0){
	print( paste("Number of iterations completed=", s) )
	print( paste("k=", dn) )
	#par(mfrow=c(1,1))
	#k_s=max(config_s)
	#plot(x[,1],x[,2])
	#for(i in 1:k_s){
	#points(x[config_s==i,1],x[config_s==i,2],col=i)
	#}
	}


	#if s is bigger than burnin, save the output
	if (s>=burnin){
		config[s+1-burnin,]=config_s
		mu_x[[s+1-burnin]]=mu_x_s
		sigma_x[[s+1-burnin]]=sigma_x_s	
	}

}

#return output
output=list(config = config, mu_x=mu_x, sigma_x=sigma_x)
return( output)

}
