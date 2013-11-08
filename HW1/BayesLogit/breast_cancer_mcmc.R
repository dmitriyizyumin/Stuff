
##
#
# Logistic regression
# 
# Y_{i} | \beta \sim \textrm{Bin}\left(n_{i},e^{x_{i}^{T}\beta}/(1+e^{x_{i}^{T}\beta})\right)
# \beta \sim N\left(\beta_{0},\Sigma_{0}\right)
#
##

library(MASS)
library(mvtnorm)
library(coda)

########################################################################################
########################################################################################
## Handle batch job arguments:

# 1-indexed version is used now.
args <- commandArgs(TRUE)

cat(paste0("Command-line arguments:\n"))
print(args)

####
# sim_start ==> Lowest simulation number to be analyzed by this particular batch job
###

#######################
sim_start <- 1000
length.datasets <- 200
#######################

if (length(args)==0){
  sinkit <- FALSE
  sim_num <- sim_start + 1
  set.seed(1330931)
} else {
  # Sink output to file?
  sinkit <- TRUE
  # Decide on the job number, usually start at 1000:
  sim_num <- sim_start + as.numeric(args[1])
  # Set a different random seed for every job number!!!
  set.seed(762*sim_num + 1330931)
}

# Simulation datasets numbered 1001-1200

########################################################################################
########################################################################################

"bayes.logreg" <- function(m,y,X,beta.0,Sigma.0.inv,niter=10,burnin=10,
                           print.every=1000,retune=400,verbose=TRUE,v)
{
	#MH to be used within Gibbs
	"metropolis" <- function(theta,pos,v.j)
	{
		#Logit Inverse Calculator used in Log Stationary
		"logitinv" = function(x,theta.c)
		{
			u = t(x) %*% theta.c
			return(exp(u)/(1+exp(u)))
		}

		#Use to Calculate Log of Stationary Distributions
		"log.stationary" = function(theta.i)
		{
			n = length(y)
			store = rep(0,n)
			add = (-1/2)*theta.i^2

			for(i in 1:n)
			{
				theta[pos] = theta.i
				li = logitinv(X[i,],theta)
				one = y[i]*log(li); two = (m[i]-y[i])*(log(1-li))
			
				if(is.nan(one)) 
					one = 0
			
				if(is.nan(two)) 
					two = 0

				store[i] = one + two
			}

			return(sum(store)+add)
		}

		#Set values of Theta.t and Theta.star
		theta.t = theta[pos]; sigma = v.j
		theta.star = rnorm(1,mean=theta.t,sd=sigma)

		#Set up Necessary Components for Alpha
		U = runif(1,min=0,max=1)
		ls.theta.star = log.stationary(theta.star)
		ls.theta.t = log.stationary(theta.t)

		#If Accepted, then change Theta.t to be Theta.star
		if(log(U) < (ls.theta.star - ls.theta.t))
		{
			theta.t = theta.star
		}

		#Return the Result (Changed or Unchanged)
		return(theta.t)
	}

	#Set up Gibbs Sampler	
	p = 10; j = 1; N = niter + burnin
	theta = rep(0,p); count = rep(0,p)
	store = matrix(0,p,N); accept = matrix(0,p,length(v))
	
	#Gibbs Sampler using Metropolis
	for(t in 1:N)
	{
		#Store the Old Parameters
		old = theta

		#Call MH on each Theta and Store Result
		for(i in 1:10)
		{
			theta[i] = metropolis(theta,i,v[i])
			store[i,t] = theta[i]
		}

		#Increment Count (for Printing)
		for(i in 1:10)
		{
			if(theta[i] != old[i])
				count[i] = count[i] + 1
		}
	}

	for(i in 1:p)
	{
		cat("Acceptance rate for v",i,"is",count[i]/N,"\n")
	}

	#store.new = store[,(burnin+1):(niter+burnin)]
	store.new = store
	return(store.new)
}



#################################################
# Set up the specifications:
beta.0 <- matrix(c(0,0))
p <- 2; Sigma.0.inv <- diag(rep(1.0,p))
niter <- 500; burnin <- 0
# etc... (more needed here)
#################################################

# Read data corresponding to appropriate sim_num:
setwd("C:/Users/Justin/Desktop/STA 250")

fileName = "breast_cancer.txt"
data = read.table(fileName,header=TRUE)

# Extract X and y:
X = as.matrix(data[1:10])
y = ifelse(data[11]=="M",1,0); m = rep(1,length(y))

v = rep(0,10)

for(i in 1:10)
{
	v[i] = sd(X[,i])
}

# Fit the Bayesian model:
result <- bayes.logreg(m,y,X,beta.0,Sigma.0.inv,niter=niter,burnin=burnin,v=v)

# Extract posterior quantiles...
result.mat = mcmc(t(result))
summary(result.mat)
plot(result.mat)

seq.1 = seq(0.01,0.99,by=0.01)
theta.1.quantiles = quantile(result[1,],probs=c(seq.1))
theta.2.quantiles = quantile(result[2,],probs=c(seq.1))

# Write results to a (99 x p) csv file...
write = cbind(theta.1.quantiles,theta.2.quantiles)
write.table(write,paste("results/blr_res_",sim_num_c,".csv",sep=""),sep=",",
row.names=FALSE,col.names=FALSE)

# Go celebrate.
 
cat("done. :)\n")
