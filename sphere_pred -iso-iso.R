rm(list = ls())
source("func.R")
set.seed(1)


######True parameters for gamma functions######
# parms: 
alpha<- c(-0.5,0,0)  
beta <- c(-0.5,0,0)  
kappa <- 0
nu=0.5
range= -2/log(0.05)
nuggets = 0.05^2 
# locs
nsrt <- 30
ns <- nsrt * nsrt
grd.obj <- grd_gen(nsrt)
prop.train = 0.7
mask.train = mask_gen_rnd(ns, prop.train)
mask.test = !mask.train
locxyz.all = cart(grd.obj$grd.all[, 1], grd.obj$grd.all[, 2])
# generate obs
z.all <- z_gen(alpha, beta, grd.obj$grd.all, kappa, nu, range, nuggets)
# divide obs
locs <- locxyz.all[mask.train, ]
locs.pred <- locxyz.all[mask.test, ]
z <- z.all[mask.train]
# Estimate parms
m.use <- 10
par0 <- c(-1, 0, 0, -1, 0, 0, kappa, nu, range)
par.mask <- rep(F, length(par0))
par.mask[c(1, 4)] <- T
n.MCMC <- 5000
burnin <- 1000
par.est <- parm_est(par0, par.mask, z, locs, m.use, nuggets, n.MCMC, burnin)
# Pred unknown locs
z.pred <- resp_pred(par.est, z, locs, locs.pred, m.use, nuggets)[[1]]
# Plot known locs
sphere_plot(z.all, rep(T, length(z.all)), grd.obj$lon, grd.obj$lat)






