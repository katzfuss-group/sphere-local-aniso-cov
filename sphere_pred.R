rm(list = ls())
library(fields)
library(GPvecchia)
library(adaptMCMC)
library(scoringRules)

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

# axially
#alpha<- c(-0.5, 0, 1.44)  
#beta <- c(-3.2, 0, 1.44) 

# iso
#alpha<- c(-0.5, 0, 0)  
#beta <- c(-0.5, 0, 0)  

alpha<- c(-0.5, 0, 1.44)  
beta <- c(-3.2, 0, 1.44)    
kappa <- 0
nu=0.5
range= -2/log(0.05)
nuggets = 0.05^2

ns <- 50*50 # grid total observation
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

set.seed(1)
z.all <- t(chol(CwithNug)) %*% rnorm(nrow(grd.all))


#########################################
locs = locxyz.all[mask.train, ]
z = z.all[mask.train]
covparms = c(alpha, beta, kappa, nu, range)
nuggets = nuggets

#########################################
## vecchia likelihood plot ##
loglikvec = c()
for(m in 1:30){
  m.use = m
  vecchia.approx=vecchia_specify(locs = locs, m = m.use, cond.yz = 'y',ic0=TRUE)
  loglikvec[m] = vecchia_likelihood(z = z, vecchia.approx = vecchia.approx, covparms = covparms, nuggets = nuggets, covmodel = "sphere")
}

plot(loglikvec ~ c(1:30), xlab="m", ylab="log-likelihood", type="o")

#########################################
## prediction ##
m = 15
vecchia.approx=vecchia_specify(locs = locs, m = m, cond.yz = 'y',ic0=TRUE)
vecchia_likelihood(z = z, vecchia.approx = vecchia.approx, covparms = covparms, nuggets = nuggets, covmodel = "sphere")

locs.pred = locxyz.all[mask.test, ]
vecchia.approx=vecchia_specify(locs = locs, locs.pred=locs.pred, m = m, cond.yz = 'y',ic0=TRUE)
preds=vecchia_prediction(z = z, vecchia.approx = vecchia.approx, covparms = covparms, nuggets = nuggets, covmodel = "sphere")
preds$mu.pred

#par(mfrow=c(1,1))
#zlim = c(min(z.all),max(z.all))
#quilt.plot(grd.all[mask.train,1], grd.all[mask.train,2],z,zlim=zlim,nx=30,ny=30)
#quilt.plot(grd.all[mask.train,1], grd.all[mask.train,2],preds$mu.obs,zlim=zlim,nx=30,ny=30)
#quilt.plot(grd.all[mask.test,1], grd.all[mask.test,2],z.all[mask.test],zlim=zlim,nx=30,ny=30)
#quilt.plot(grd.all[mask.test,1], grd.all[mask.test,2],preds$mu.pred,zlim=zlim,nx=30,ny=30)


################ Parameter estimation MCMC ##############################################
set.seed(1)
m.use=10  #conditioning set size
vecchia.approx = vecchia_specify(locs = locs, m = m.use, cond.yz = 'y',ic0=TRUE)

####### estimate All parameter ################################
######## params = c(alpha, beta, kappa, nu, range, nuggets)
####### loglikMCMC(c(alpha, beta, kappa, nu, range, nuggets))
# loglikMCMC = function(par){
#   par.pos = par
#   par.pos[8:10]=abs(par.pos[8:10])
#   return(vecchia_likelihood(
#     z = z.all[mask.train],
#     vecchia.approx = vecchia.approx,
#     # alpha, beta, kappa, nu, range
#     covparms = par.pos[1:9],
#     nuggets = par.pos[10],
#     covmodel = "sphere"
#   ))
# }

###### estimate slope only #################################
loglikMCMC = function(par){
  return(vecchia_likelihood(
    z = z.all[mask.train],
    vecchia.approx = vecchia.approx,
    # alpha, beta, kappa, nu, range
    covparms = c(par, kappa,nu,range),
    nuggets = nuggets,
    covmodel = "sphere"
  ))
}

n.MCMC = 5000

samp = MCMC(p=loglikMCMC, n=n.MCMC, init=c(-1, 0.5, 0.5, -1, 0.5, 0.5), scale = rep(1, 6),
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
par("mar")
par(mar=c(1,1,1,1))
plot(samp.coda)

## Autocorrelation function for Markov chains
autocorr(samp.coda)
#autocorr.plot(samp.coda)
autocorr.diag(samp.coda)

## Cross correlations for MCMC output
crosscorr(samp.coda)
crosscorr.plot(samp.coda,  col = topo.colors(15))

## the effective number of independent samples 
effectiveSize(samp.coda)
##################################################################

burnin=1000
#plot(acf(samp$samples[burnin:n.MCMC,]))
samp.eff=samp$samples[seq(burnin, n.MCMC,by=10),] #thinned values after convergence
par.est = colMeans(samp.eff)
par.est





