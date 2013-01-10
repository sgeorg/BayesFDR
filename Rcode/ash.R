#return ABF for vector of betahat and standard errors
ABF = function(betahat, sebetahat,sigmaa){
T = betahat/sebetahat
lambda = sebetahat^2/(sebetahat^2+sigmaa^2)
return((sqrt(lambda) * exp(0.5*T^2 *(1-lambda))))
}

logABF = function(betahat,sebetahat,sigmaa){
T = betahat/sebetahat
lambda = sebetahat^2/(sebetahat^2+sigmaa^2)
return(0.5*log(lambda) + 0.5*T^2 *(1-lambda))
}

#return matrix of ABFs for vector of sigma-a values
#normalized by maximum of each column
#betahat is n vector, sebetahat is n vector, sigmaavec is k vector
#return is n by k matrix of ABFs
matrixABF = function(betahat, sebetahat, sigmaavec){
k = length(sigmaavec)
n = length(betahat)
labf = matrix(0,nrow=n, ncol=k)
for(i in 1:k){
labf[,i] = logABF(betahat,sebetahat,sigmaavec[i])
}
maxlabf = apply(labf, 1, max)
labf = labf - maxlabf
return(exp(labf))
}

#estimate mixture proportions of sigmaa by EM algorithm
EMest = function(betahat,sebetahat,sigmaavec,niter,pi){
abf = matrixABF(betahat,sebetahat,sigmaavec)
for(i in 1:niter){
m  = t(pi * t(abf)) 
m.rowsum = rowSums(m)
classprob = m/m.rowsum
pi = apply(classprob,2, mean)
}
return(list(pi=pi,classprob=classprob))
}

normalize = function(x){return(x/sum(x))}

#return the posterior on beta given a prior
#that is a mixture of normals (pi0,mu0,sigma0)
#and observation betahat \sim N(beta,sebetahat)
#current ABF is only for mu0=0, so would need to
#generalize that for general application
#returns list whose components are k by n matrices
#k is number of mixture components, n is number of observations
posterior_dist = function(pi0,mu0,sigma0,betahat,sebetahat){
  k= length(pi0)
  n= length(betahat)
  
  pi1 = pi0 * t(matrixABF(betahat,sebetahat,sigma0))
  pi1 = apply(pi1, 2, normalize) #pi1 is now an k by n matrix

  #make k by n matrix versions of sigma0^2 and sebetahat^2
  # and mu0 and betahat
  s0m2 = matrix(sigma0^2,nrow=k,ncol=n,byrow=F)
  sebm2 = matrix(sebetahat^2,nrow=k,ncol=n, byrow=T)
  mu0m = matrix(mu0,nrow=k,ncol=n,byrow=F)
  bhatm = matrix(betahat,nrow=k,ncol=n,byrow=T)

  sigma1 = (1/s0m2 + 1/sebm2)^(-0.5)  
  w = (1/s0m2)/(1/s0m2 + 1/sebm2)
  mu1 = w*mu0m + (1-w)*bhatm
  
  return(list(pi1=pi1,mu1=mu1,sigma1=sigma1))
}

#helper function for posterior_sample
#samples nsamp integers from 1:k according to a given prob vector
sample_component=function(p,nsamp){
  
  return(sample(length(p),nsamp,replace=T,prob=p))
}

#m is a k by n matrix
#comp is a n vector of values in 1-k
#returns the comp[i]-th row of m[,i] 
extract_component=function(comp,m){
  return(m[cbind(comp,seq(comp))])
}

#returns matrix of nsamp samples from posterior
#computed using posterior_dist
# NOTE THIS IS UNTESTED, AND PROBABLY NOT WORKING YET...
posterior_sample = function(post,nsamp){
  component = as.vector(apply(post$pi1,2,sample_component,nsamp=nsamp))
  k = ncol(post$pi1)
  s = rep(1:k,rep(nsamp,k))
  index = cbind(component,s) #set up indices of mu and sigma to extract
  m = post$mu1[index]
  ss = post$sigma1[index]
  res = matrix(rnorm(length(m),mean=m,sd=ss),nrow=nsamp)
  return(res)
}


#find point estimates of beta from a posterior produced by posterior_dist
posterior_mean = function(post){
  return(apply(post$pi1 * post$mu1,2,sum))
}

#return posterior of being >T for a mixture of Gaussians
# each of pi1, mu1, sigma1 is a k by n matrix
# jth column provides parameters for jth mixture of gauusians 
# return an n vector of probabilities
PosteriorProbExceedsT = function(pi1,mu1,sigma1,T=0){
return(apply(pi1 * pnorm(T,mu1,sigma1,lower.tail=F),2,sum))
}
  
#main adaptive shrinkage function
#takes a vector of betahats and ses;
#fits a mixture of normals to it
# and 
ash = function(betahat,sebetahat,nsamp=0){
sigmaavec = c(0.001,0.002,0.004,0.008,0.016,0.032,0.064,0.128,0.256,0.512,1.024,2.048,4.096,8.192)
pi=rgamma(length(sigmaavec),1,1)
completeobs = !is.na(betahat) & !is.na(sebetahat)
pi=EMest(betahat[completeobs],sebetahat[completeobs],sigmaavec,1000,pi)

post = posterior_dist(pi$pi,0,sigmaavec,betahat,sebetahat)
PositiveProb = PosteriorProbExceedsT(post$pi1,post$mu1,post$sigma1,0)
PosteriorMean = posterior_mean(post)
if(nsamp>0){
  sample = posterior_sample(post,nsamp)
}
return(list(post=post,sample=sample,PosteriorMean = PosteriorMean,PositiveProb =PositiveProb))
}

