library(GPvecchia)


#######   simple example from vignette

spatial.dim=2
n=50
locs <- cbind(runif(n),runif(n))

beta=2
sig2=1; range=.1; smooth=1.5
covparms =c(sig2,range,smooth)
covfun <- function(locs) sig2*MaternFun(fields::rdist(locs),covparms)
nuggets=rep(.1,n)

Om0 <- covfun(locs)+diag(nuggets)
z=as.numeric(t(chol(Om0))%*%rnorm(n))

n.p=100
grid.oneside=seq(0,1,length=round(sqrt(n.p)))
locs.pred=as.matrix(expand.grid(grid.oneside,grid.oneside)) # grid of pred.locs
n.p=nrow(locs.pred)


m=30
##  specify predictions --- edited from vignette to use cond.yz='y'
vecchia.approx=vecchia_specify(locs,m,locs.pred=locs.pred,cond.yz='y',ic0=TRUE)
preds=vecchia_prediction(z,vecchia.approx,covparms,nuggets)





#######  drawing samples from posterior predictive distribution (PP)

# desired number of draws from PP
N=1000 

# index vector to extract and order prediction locations
orig.order=order(preds$U.obj$ord)
ord.pred=seq(n+n.p,1,by=-1)[orig.order[!preds$U.obj$obs[orig.order]]]

# draws from standard normal
snd=matrix(rnorm(N*(n+n.p)),ncol=N)

# transform draws to PP 
post.draws=Matrix::solve(Matrix::t(preds$V.ord),snd)[ord.pred,]+preds$mu.pred



#### check that draws roughly agree with exact PP mean and var
mu.samp=apply(post.draws,1,mean)
plot(preds$mu.pred,mu.samp)

var.samp=apply(post.draws,1,var)
plot(preds$var.pred,var.samp)


