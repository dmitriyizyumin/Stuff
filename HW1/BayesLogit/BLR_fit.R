#Dmitriy Izyumin
#STA 250
#FALL 2013

##
#
# Logistic regression
# 
# Y_{i} | \beta \sim \textrm{Bin}\left(n_{i},e^{x_{i}^{T}\beta}/(1+e^{x_{i}^{T}\beta})\right)
# \beta \sim N\left(\beta_{0},\Sigma_{0}\right)
#
##

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

#mydata<-read.csv('C:/Users/Dmitriy/Documents/GitHub/Stuff/HW1/BayesLogit/data/blr_data_1111.csv',header=T)


"bayes.logreg" <- function(n,y,X,beta.0,Sigma.0.inv,niter=10000,burnin=1000,
                           print.every=1000,retune=100,verbose=TRUE)
{
  
  sigma<-diag(solve(Sigma.0.inv))
  beta.t<-beta.0
  
  total.iter<-burnin+niter
  retune.times<-numeric(total.iter)
  retune.times[seq(from=retune,by=retune,length.out=floor(burnin/retune))]<-1
  accept.nums<-numeric(2)
  draws<-matrix(NA,nrow=niter,ncol=2)
  
  for(t in 1:total.iter){
    
    for(i in 1:2){  
      beta.star<-beta.t
      beta.star[i]<-rnorm(1,beta.t[i],sigma[i])
      log.pi.dif<-(t(y)%*%X%*%beta.star-n%*%log(1+exp(X%*%beta.star)))-(t(y)%*%X%*%beta.t-n%*%log(1+exp(X%*%beta.t)))
      log.alpha<-min(0,log.pi.dif)
      
      #To catch NAs
      if(is.na(log.alpha)){
        log.alpha <- -Inf
      }
      
      #Accept proposal with probability alpha                                                
      if(log(runif(1,0,1))<log.alpha){
        beta.t<-beta.star
        accept.nums[i]<-accept.nums[i]+1
      }
                                                           
    } #end iteration
    
    #record draws
    if(t>burnin){
      draws[t-burnin,]<-beta.t
    }
    
    #retune sigma
    if(retune.times[t]){
      accept.rates<-accept.num/retune
      sigma<-sigma*(1 + 0.3*(accept.rates<0.3) - 0.3*(accept.rates>0.6))
      accept.nums<-numeric(2)      
    }
  
  } #end process
  return(draws)
  
} #end function bayes.logreg

#################################################
# Set up the specifications:
beta.0 <- matrix(c(0,0))
Sigma.0.inv <- diag(rep(1.0,2))
niter <- 10000
burnin=5000
# etc... (more needed here)
#################################################

# Read data corresponding to appropriate sim_num:

# Extract X and y:
data<-read.table(file=paste0('data/blr_data_',as.character(sim_num),'.csv'),sep=',',header=T)

n<-data[,2]
X<-data[,3:4]
y<-data[,1]

# Fit the Bayesian model:
mydraws<-bayes.logreg(n,y,X,beta.0,Sigma.0.inv,verbose=T)

# Extract posterior quantiles...
q.beta<-cbind(quantile(mydraws[,1],probs=seq(0.01,0.99,0.01)),quantile(mydraws[,2],probs=seq(0.01,0.99,0.01)))


# Write results to a (99 x p) csv file...
write.table(x=q.beta,file=paste0('results/blr_res_',as.character(sim_num),'.csv'),sep=',',col.names=F,row.names=F)

# Go celebrate.
 
cat("done. :)\n")







