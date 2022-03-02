rm(list = ls())
library(fields)
library(GPvecchia)
library(adaptMCMC)
library(scoringRules)

set.seed(1)

###### Define functions for rotation matrices ######
Rx <- function(x){
  return(matrix(c(1,0,0,0,cos(x),-sin(x),0,sin(x),cos(x)),3,3,byrow=TRUE))
}

Ry <- function(y){
  return(matrix(c(cos(y),0,sin(y),0,1,0,-sin(y),0,cos(y)),3,3,byrow=TRUE))
}

Rz <- function(z){
  return(matrix(c(cos(z),-sin(z),0,sin(z),cos(z),0,0,0,1),3,3,byrow=TRUE))
}

######Define the function for transform (l,L)-coordinates into (x,y,z)-coordinates######
###### (l=longitude, L=Latitude)
cart <- function(l,L){  
  x = cos(L)*cos(l)
  y = cos(L)*sin(l)
  z = sin(L)
  return(c(x,y,z))
}


################################################################################################
################################################################################################
######True parameters for gamma functions######

## sin version: 
alpha<- c(-0.5, 0, 0)  
beta <- c(-0.5, 0, 0) 
kappa <- 0
nu=0.5
range= -2/log(0.05)
nuggets = 0.05^2

ns <- 30*30 # grid total observation
lon=seq(-pi+0.1, pi-0.1, length.out = sqrt(ns))
lat=seq(-pi/2+0.05, pi/2-0.05, length.out = sqrt(ns))

#grd <- cbind(runif(ns,min=-pi,max=pi),runif(ns,min=-pi/2,max=pi/2))
grd.all = as.matrix(expand.grid(lon,lat))

prop.train = 0.7
index.train = sample(1:ns, size =prop.train*ns, replace =FALSE)
mask.train=rep(FALSE, ns)
mask.train[index.train]=TRUE
mask.test=!mask.train

locxyz.all = matrix(NA, nrow(grd.all), 3)
for(i in 1:nrow(grd.all)) {
  locxyz.all[i,] = cart(grd.all[i,1],grd.all[i,2])
}

# sin version
gamma1 <- exp(alpha[1]+alpha[2]*sin(grd.all[,1])+alpha[3]*grd.all[,2])
gamma2 <- exp(beta[1]+beta[2]*sin(grd.all[,1])+beta[3]*grd.all[,2])
Sigma <- list()

for (i in 1:nrow(grd.all)){
  D <- diag(c(1,gamma1[i],gamma2[i]))
  Sigmatilde <- Rx(kappa)%*%D%*%t(Rx(kappa))
  Sigma[[i]] <- Rz(grd.all[i,1])%*%Ry(grd.all[i,2])%*%Sigmatilde%*%t(Ry(grd.all[i,2]))%*%t(Rz(grd.all[i,1]))
}

######Calculate Mahalanobis distance and normalization term c(,)
dis.mat <- c.mat <- matrix(0,nrow(grd.all),nrow(grd.all))
for (i in 1:nrow(grd.all)){
  for (j in 1:i){
    dis.mat[i,j] <- (2*t(cart(grd.all[i,1],grd.all[i,2])-cart(grd.all[j,1],grd.all[j,2]))%*%solve(Sigma[[i]]+Sigma[[j]])%*%(cart(grd.all[i,1],grd.all[i,2])-cart(grd.all[j,1],grd.all[j,2])))^(1/2)
    dis.mat[j,i]=dis.mat[i,j]
    c.mat[i,j] <- (det(Sigma[[i]]))^(1/4) * (det(Sigma[[j]]))^(1/4) * (det((Sigma[[i]]+Sigma[[j]])/2))^(-1/2)
    c.mat[j,i]=c.mat[i,j] 
  }
  print(i)
}

C <- c.mat * Matern(dis.mat,range=range,nu=nu) 
CwithNug = C + diag(nuggets,nrow=nrow(grd.all),ncol=nrow(grd.all))

z.all <- t(chol(CwithNug)) %*% rnorm(nrow(grd.all))


#########################################
locs = locxyz.all[mask.train, ]
locs.pred = locxyz.all[mask.test, ]
z = z.all[mask.train]
covparms = c(alpha, beta, kappa, nu, range)
nuggets = nuggets


################ Parameter estimation MCMC ##############################################
m.use=10  #conditioning set size
vecchia.approx = vecchia_specify(locs = locs, m = m.use, cond.yz = 'y',ic0=TRUE)


###### estimate slope only, use axially symmetric model #################################
## alpha =c( par1, 0,par2)
## beta = c( par3,0,par4)
## par = c(par1, par2,par3,par4)
loglikMCMC = function(par){
  return(vecchia_likelihood(
    z = z.all[mask.train],
    vecchia.approx = vecchia.approx,
    # alpha, beta, kappa, nu, range
    covparms = c(par[1],0,par[2],par[3],0,par[4], kappa,nu,range),
    nuggets = nuggets,
    covmodel = "sphere"
  ))
}

n.MCMC = 5000

samp = MCMC(p=loglikMCMC, n=n.MCMC, init=c(-1,0.2, -1,0.2), scale = rep(1, 4),
            adapt = TRUE, acc.rate = 0.23, gamma = 2/3,
            showProgressBar=TRUE)

samp0=samp
str(samp)
summary(samp$samples)
###########################################################################
## covariance of last jump distribution
samp$cov.jump

## convert chans into 'coda' object
samp.coda <- convert.to.coda(samp)

## plot chains and marginals
require(coda)
#par("mar")
#par(mar=c(1,1,1,1))
#plot(samp.coda)

## Autocorrelation function for Markov chains
#autocorr(samp.coda)
#autocorr.plot(samp.coda)
#autocorr.diag(samp.coda)

## Cross correlations for MCMC output
#crosscorr(samp.coda)
#crosscorr.plot(samp.coda,  col = topo.colors(15))

## the effective number of independent samples 
#effectiveSize(samp.coda)
##################################################################

burnin=1000
samp.eff=samp$samples[seq(burnin, n.MCMC,by=10),] #thinned values after convergence
par.est = colMeans(samp.eff)
par.est   






