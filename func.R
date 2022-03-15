library(parallel)
library(fields)
library(GPvecchia)
library(adaptMCMC)
library(scoringRules)
library(oce)


###### Define functions for rotation matrices ######
Rx <- function(x) {
    return(matrix(c(
        1, 0, 0, 0, cos(x), -sin(x), 0, sin(x), cos(x)
    ), 3, 3, byrow = TRUE))
}

Ry <- function(y) {
    return(matrix(c(
        cos(y), 0, sin(y), 0, 1, 0, -sin(y), 0, cos(y)
    ), 3, 3, byrow = TRUE))
}

Rz <- function(z) {
    return(matrix(c(
        cos(z), -sin(z), 0, sin(z), cos(z), 0, 0, 0, 1
    ), 3, 3, byrow = TRUE))
}


######Define the function for transform (l,L)-coordinates into (x,y,z)-coordinates######
###### (l=longitude, L=Latitude)
cart <- function(l, L) {
    x = cos(L) * cos(l)
    y = cos(L) * sin(l)
    z = sin(L)
    return(cbind(x, y, z))
}


###### Generate lon/lat grid ######
grd_gen <- function(nsrt) {
    ns <- nsrt * nsrt # grid total observation
    lon = seq(-pi + 0.1, pi - 0.1, length.out = nsrt)
    lat = seq(-pi / 2 + 0.05, pi / 2 - 0.05, length.out = nsrt)
    
    grd.all = as.matrix(expand.grid(lon, lat))
    return(list(lon = lon, lat = lat, grd.all = grd.all))
}


###### Generate mask ######
mask_gen <- function(ns, prop) {
    index = sample(1 : ns, size = prop * ns, replace = FALSE)
    mask = rep(FALSE, ns)
    mask[index] = TRUE
    return(mask)
}


###### Generate responses ######
z_gen <- function(gamma1, gamma2, grd.all, kappa, nu, range, nuggets, 
                  cluster = NULL){
    n <- nrow(grd.all)
    if(is.null(cluster)){
        inner.cluster = T
        cluster <- makeCluster(detectCores())
    }

    Sigma.input <- cbind(grd.all, gamma1, gamma2, kappa)
    Sigma.input <- lapply(1 : nrow(Sigma.input), 
                          function(x){return(Sigma.input[x,])})
    Sigma.func <- function(x){
        D <- diag(c(1,x[3],x[4]))
        Sigmatilde <- Rx(x[5]) %*% D %*% t(Rx(x[5]))
        Rz(x[1]) %*% Ry(x[2]) %*% Sigmatilde %*% t(Ry(x[2])) %*% 
            t(Rz(x[1]))
    }
    clusterExport(cluster, c("Rx", "Ry", "Rz", "cart"), envir = .GlobalEnv)
    Sigma <- parLapply(cluster, Sigma.input, Sigma.func) 
    clusterExport(cluster, c("grd.all", "Sigma", "n"), envir = environment())
    dis.func <- function(x){
        if(x < 1)
            return(0)
        i = (x - 1) %% n + 1
        j = floor((x - 1) / n) + 1
        cart.diff = cart(grd.all[i,1], grd.all[i,2])-
            cart(grd.all[j,1], grd.all[j,2])
        2 * (cart.diff %*% solve(Sigma[[i]] + Sigma[[j]]) %*% 
                 t(cart.diff))^(1/2)
    }
    c.func <- function(x){
        if(x < 1)
            return(0)
        i = (x - 1) %% n + 1
        j = floor((x - 1) / n) + 1
        (det(Sigma[[i]]))^(1/4) * (det(Sigma[[j]]))^(1/4) * 
            (det((Sigma[[i]]+Sigma[[j]])/2))^(-1/2)
    }
    tmp.mat <- matrix(1 : n, n, n, F) + outer(rep(n, n), 0 : (n - 1))
    tmp.mat[lower.tri(tmp.mat)] <- 0
    dis.mat <- parApply(cl = cluster, X = tmp.mat, MARGIN = c(1, 2),
                        FUN = dis.func)
    # dis.mat <- apply(X = tmp.mat, MARGIN = c(1, 2), FUN = dis.func)
    c.mat <- parApply(cl = cluster, X = tmp.mat, MARGIN = c(1, 2), 
                        FUN = c.func)
    C <- c.mat * Matern(dis.mat, range=range, nu=nu) 
    CwithNug = C + diag(nuggets, nrow=n, ncol=n)
    z.all <- t(chol(CwithNug)) %*% rnorm(n)
    
    if(inner.cluster)
        stopCluster(cluster)
    
    z.all
}


###### Parameter Estimation ######
parm_est <- function(par0, mask, z, locs, m, nuggets, n.MCMC, burnin)
{
    vecchia.approx <- vecchia_specify(locs = locs, m = m, 
                                      cond.yz = 'y',ic0=TRUE)
    loglikMCMC <- function(par){
        parLong <- par0
        parLong[mask] <- par
        return(vecchia_likelihood(
            z = z,
            vecchia.approx = vecchia.approx,
            # alpha, beta, kappa, nu, range
            covparms = parLong,
            nuggets = nuggets,
            covmodel = "sphere"
        ))
    }
    scale <- rep(1, sum(mask))
    samp <- MCMC(p=loglikMCMC, n=n.MCMC, init=par0[mask], scale = scale, 
                 adapt = TRUE, acc.rate = 0.23, gamma = 2/3,
                 showProgressBar=TRUE)
    #thinned values after convergence
    samp.eff <- samp$samples[seq(burnin, n.MCMC, by=10), ] 
    #par.est
    colMeans(samp.eff)
}


###### Prediction Unknown Locs ######
resp_pred <- function(par, z, locs, locs.pred, m, nuggets)
{
    vecchia.approx <- vecchia_specify(locs = locs, locs.pred=locs.pred, m = m, 
                                      cond.yz = 'y',ic0=TRUE)
    preds <- vecchia_prediction(z = z, vecchia.approx = vecchia.approx, 
                                covparms = par, nuggets = nuggets, 
                                covmodel = "sphere")
    preds$mu.pred
}


###### Prediction Unknown Locs ######
sphere_plot <- function(z.all, mask, lon, lat, zlim = range(z.all, na.rm = T)){
    z.all[!mask] <- NA
    z.all.mat <- matrix(z.all, nrow = length(lon))
    cm <- colormap(zlim = zlim)
    par(mar = c(2, 1, 1, 1))
    drawPalette(colormap = cm)
    mapPlot(longitude = c(-180, 180, 0, 0), latitude = c(0, 0, -90, 90), 
            grid = TRUE, col = "lightgray", drawBox = FALSE, 
            longitudelim = c(-180, 180), type = "n") # defaults to moll projection
    mapImage(longitude = lon*180/pi, latitude = lat*180/pi, z = z.all.mat, 
             zlim = zlim, colormap = cm, missingColor = NA)
}






