#Dmitriy Izyumin
#Sta 250
#Fall 2013

#HW1 Q3

cancer.data<-read.table('breast_cancer.txt',header=T)

n<-dim(cancer.data)[1]
y<-cancer.data[,11]
y<-as.numeric(y=='M')
X<-as.matrix(cancer.data[,1:10])
X<-cbind(rep(1,n),scale(X))
p<-dim(X)[2]

#initial values
beta.0<-rep(0,p)
sigma.0<-rep(1,p)

#specifics of the MC
burnin<-5000
niter<-10000
retune<-500
retune.times<-numeric(burnin+niter)
retune.times[seq(from=retune,by=retune,length.out=floor(burnin/retune))]<-1

beta.t<-beta.0
sigma<-diag(diag(sigma.0))
n.accept<-numeric(p)
draws<-matrix(NA,nrow=niter+burnin,ncol=p)

for(t in 1:(niter+burnin)){
  
  for(i in 1:p){  
    beta.star<-beta.t
    beta.star[i]<-rnorm(1,beta.t[i],sigma[i])
    log.pi.dif<-(t(y)%*%X%*%beta.star-sum(log(1+exp(X%*%beta.star))))-(t(y)%*%X%*%beta.t-sum(log(1+exp(X%*%beta.t))))
    log.alpha<-min(0,log.pi.dif)
    
    #To catch NAs
    if(is.na(log.alpha)){
      log.alpha <- -Inf
    }
    
    #Accept proposal with probability alpha                                                
    if(log(runif(1))<log.alpha){
      beta.t<-beta.star
      n.accept[i]<-n.accept[i]+1
    }                                                          
  } #end iteration
  
  #record draws
  draws[t,]<-beta.t
  
  #retune sigma
  if(retune.times[t]){
    accept.rates<-n.accept/retune
    cat(paste0('STEPS ',t-retune+1,' THROUGH ',t,':\n'))
    for(i in 1:p){
      cat(paste0('Acceptance rate for beta',i,' was ',100*round(accept.rates[i],2),"%\n"))     
    }
    sigma<-sigma * (1 - 0.4*(accept.rates<0.3) + 0.4*(accept.rates>0.6))
    n.accept<-numeric(p)      
  }
  
}#end process

#Delete burnin
draws<-draws[-(1:burnin),]


#lag-1 autocorrelations
sapply(1:11,function(i){acf(draws[i,],plot=F)[1]$acf})


#posterior predictive check
pred.samps<-matrix(NA,nrow=niter,ncol=n)
for(i in 1:niter){
  theta<-draws[1,]
  pred.samps[i,]<-sapply(X=exp(X%*%theta)/(1+exp(X%*%theta)),FUN=function(p){rbinom(1,1,p)})  
}
means<-rowMeans(pred.samps)
par(mfrow=c(1,1))
plot(density(means),main='Distribution of the Means of Predictive Samples')
abline(v=mean(y),col='red')
2*min(sum(means<mean(y)),sum(means>mean(y)))/niter  #p-value
