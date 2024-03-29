library(parallel)
library(fields)
library(GPvecchia)
library(adaptMCMC)
library(scoringRules)
library(oce)
library(ncdf4)


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
mask_gen_rnd <- function(ns, prop) {
    index <- sample(1 : ns, size = prop * ns, replace = FALSE)
    mask <- rep(FALSE, ns)
    mask[index] <- TRUE
    return(mask)
}
mask_gen_region <- function(grd.all, lon.width, lat.width) {
    ns <- nrow(grd.all)
    index <- sample(1 : ns, size = 1)
    lon <- grd.all[index, 1]
    lat <- grd.all[index, 2]
    
    lon.mask <- abs(grd.all[, 1] - lon) < lon.width | 
        abs(grd.all[, 1] - lon - 2*pi) < lon.width |
        abs(grd.all[, 1] - lon + 2*pi) < lon.width
    lat.mask <- abs(grd.all[, 2] - lat) < lat.width 
    return(!(lon.mask & lat.mask))
}
mask_gen_Nile <- function(grd.all) {
  lon <- grd.all[, 1]
  lat <- grd.all[, 2]
  lon.mask <- (lon >= 29.9 / 180 * pi) & (lon <= 32.6 / 180 * pi)
  lat.mask <- (lat >= 22 / 180 * pi) & (lat <= 30 / 180 * pi)
  return(!(lon.mask & lat.mask))
}

###### Generate responses ######
z_gen <- function(alpha, beta, grd.all, kappa, nu, range, nuggets, 
                  cluster = NULL){
    n <- nrow(grd.all)
    gamma1 <- exp(alpha[1]+alpha[2]*sin(grd.all[,1])+alpha[3]*grd.all[,2])
    gamma2 <- exp(beta[1]+beta[2]*sin(grd.all[,1])+beta[3]*grd.all[,2])
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
        sqrt(2) * (cart.diff %*% solve(Sigma[[i]] + Sigma[[j]]) %*% 
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
parm_est <- function(par0, mask, z, locs, m, nuggets, n.MCMC, burnin,
                     debugFn = NULL)
{
    vecchia.approx <- vecchia_specify(locs = locs, m = m,
                                      cond.yz = 'y',ic0=TRUE)
    loglikMCMC <- function(par){
        parLong <- par0
        parLong[mask] <- par
        # special treatment for range
        # the last coef in par0 is treated as log range if 
        #   the last coef in mask if TRUE
        if(mask[length(mask)])
            parLong[length(mask)] <- exp(parLong[length(mask)])
        # special treatment for iso
        #   set beta10 and beta20 to be the same
        if(mask[1] == T & mask[4] == F)
            parLong[4] <- parLong[1]
        cat("Parm = ", parLong)
        if(any(abs(parLong) > 20)){ # reduce this value if singularity occurs
            loglk <- -Inf
        }else{
            # Vecchia loglk
            loglk <- tryCatch(vecchia_likelihood(
                z = z,
                vecchia.approx = vecchia.approx,
                # alpha, beta, kappa, nu, range
                covparms = parLong,
                nuggets = nuggets,
                covmodel = "sphere"
            ), warning = function(x){return(-Inf)},
            error = function(x){return(-Inf)})
        }
        if(is.nan(loglk))
            loglk <- -Inf
        cat(" loglk = ", loglk, "\n")
        return(loglk)
    }
    scale <- rep(1, sum(mask))
    samp <- MCMC(p=loglikMCMC, n=n.MCMC, init=par0[mask], scale = scale,
                 adapt = TRUE, acc.rate = 0.23, gamma = 2/3,
                 showProgressBar=TRUE)
    # debug statement
    if(!is.null(debugFn))
        save(samp, file = debugFn)
    #thinned values after convergence
    samp.eff <- samp$samples[seq(burnin, n.MCMC, by=10), , drop = F]
    #par.est
    par.est <- matrix(par0, nrow(samp.eff), length(par0), byrow = T)
    par.est[, mask] <- samp.eff
    par.est
}


###### Prediction Unknown Locs ######
resp_pred <- function(par, z, locs, locs.pred, m, nuggets)
{
    vecchia.approx <- vecchia_specify(locs = locs, locs.pred=locs.pred, m = m, 
                                      cond.yz = 'y',ic0=TRUE)
    preds <- vecchia_prediction(z = z, vecchia.approx = vecchia.approx, 
                                covparms = par, nuggets = nuggets, 
                                covmodel = "sphere")
    return(preds)
}


###### Plot sphere projected to 2D ######
sphere_plot <- function(z.all, mask, lon, lat, zlim = range(z.all, na.rm = T), 
                        draw.palette = T, fn = NULL, fig.width = 7, 
                        fig.height = 5){
    z.all[!mask] <- NA
    z.all.mat <- matrix(z.all, nrow = length(lon))
    cm <- colormap(zlim = zlim)
    par(mar = c(2, 1, 1, 1))
    if(!is.null(fn))
        pdf(fn, fig.width, fig.height)
    if(draw.palette)
        drawPalette(colormap = cm, cex = 2.5)
    mapPlot(longitude = c(-180, 180, 0, 0), latitude = c(0, 0, -90, 90), 
            grid = TRUE, col = "lightgray", drawBox = FALSE, 
            longitudelim = c(-180, 180), type = "n") # defaults to moll projection
    mapImage(longitude = lon*180/pi, latitude = lat*180/pi, z = z.all.mat, 
             zlim = zlim, colormap = cm, missingColor = NA)
    if(!is.null(fn))
        dev.off()
}


###### Compute MAE, RMSE, CRPS, Energy ######
score_func <- function(z, preds, n.energy = 1000){
    MAE <- mean(abs(z - colMeans(preds)))
    RMSE <- sqrt(mean((z - colMeans(preds))^2))
    CRPS <- mean(crps_norm(z, mean = colMeans(preds), sd = apply(preds, 2, sd)))
    Energy <- es_sample(z, t(preds))
    ret <- c(MAE, RMSE, CRPS, Energy)
    names(ret) <- c("MAE", "RMSE", "CRPS", "Energy")
    return(ret)
}


###### Function for sim/real application studies  ######
sim_func <- function(ns, grd.obj, z.all, nu, range, nuggets, 
                     type, plot.flag = F, seed = 1){
    prop.train <- 0.8
    n.test.region <- 10
    mask.train.rnd <- mask_gen_rnd(ns, prop.train)
    mask.test.rnd <- !mask.train.rnd
    if(type == "real"){
      mask.train.region <- mask_gen_Nile(grd.obj$grd.all)
    }else{
      mask.train.region <- rep(T, ns)
      for(i in 1 : n.test.region){
        mask.train.region <- mask.train.region & 
          mask_gen_region(grd.obj$grd.all, 0.4, 0.2)
      }
    }
    mask.test.region <- !mask.train.region
    locxyz.all <- cart(grd.obj$grd.all[, 1], grd.obj$grd.all[, 2])
    
    if(type == "real"){
        par0.lst <- list(c(-0.2736117,.0,.0,-1.5002264,.0,
                           .0,.0,2.500000,exp(-2.5987020)),
                         c(-0.2736117,.0,.0,-1.5002264,.0,
                           .0,.0,2.500000,exp(-2.5987020)),
                         c(-0.2736117,.0,.0,-1.5002264,.0,
                           .0,.0,2.500000,exp(-2.5987020)))
        par.mask.lst <- list(c(T, F, F, F, F, F, F, F, F),
                             c(T, F, T, T, F, T, F, F, F),
                             c(T, T, T, T, T, T, T, F, F))
    }else{
        par0.lst <- list(c(rep(0, 7), nu, range),
                         c(rep(0, 7), nu, range),
                         c(rep(0, 7), nu, range))
        par.mask.lst <- list(c(T, F, F, F, F, F, F, F, F),
                             c(T, F, T, T, F, T, F, F, F),
                             c(T, T, T, T, T, T, T, F, F))
    }
    mdl.lst <- list("iso", "axially", "nonsta")
    m <- 10
    n.MCMC <- 5000
    burnin <- 1000
    
    if(plot.flag){
        sphere_plot(z.all, mask.train.rnd, grd.obj$lon, grd.obj$lat,
                    draw.palette = F, fn = paste0("rand-", type, "-true.pdf"))
        sphere_plot(z.all, mask.train.region, grd.obj$lon, grd.obj$lat,
                    draw.palette = F, fn = paste0("rect-", type, "-true.pdf"))
    }else{
        input <- list(z.all = z.all, mask = mask.train.rnd, lon = grd.obj$lon,
                      lat = grd.obj$lat, draw.palette = F, 
                      fn = paste0("rand-", type, "-true.pdf"))
        save(input, file = paste0("rand-", type, "-true.RData"))
        input <- list(z.all = z.all, mask = mask.train.region, lon = grd.obj$lon,
                      lat = grd.obj$lat, draw.palette = F, 
                      fn = paste0("rect-", type, "-true.pdf"))
        save(input, file = paste0("rect-", type, "-true.RData"))
    }
    
    tbl <- matrix(NA, 3, 8)
    colnames(tbl) <- c(paste0(c("MAE", "RMSE", "CRPS", "Energy"), "-rnd"), 
                       paste0(c("MAE", "RMSE", "CRPS", "Energy"), "-region"))
    rownames(tbl) <- c("Isotropic", "Axially symmetric", "Nonstationary")
    tbl <- as.data.frame(tbl)
    
    for(i in 1 : length(par0.lst)){
        par0 <- par0.lst[[i]]
        par.mask <- par.mask.lst[[i]]
        # special treatment for the range par
        # change the last coef in par0 to log par0 if the last coef in
        #   par.mask is TRUE
        if(par.mask[length(par.mask)])
            par0[length(par.mask)] <- log(par0[length(par.mask)])
        parm.est.rnd <- parm_est(par0, par.mask, z.all[mask.train.rnd],
                                  locxyz.all[mask.train.rnd, ], m, nuggets,
                                  n.MCMC, burnin)
        parm.est.region <- parm_est(par0, par.mask, z.all[mask.train.region],
                                    locxyz.all[mask.train.region, ], m, nuggets,
                                    n.MCMC, burnin)
        # special treatment for the range par
        # change the last coef in parm.est to exp parm.est if the last coef in
        #   par.mask is TRUE
        if(par.mask[length(par.mask)]){
            parm.est.rnd[, length(par.mask)] <- 
                exp(parm.est.rnd[, length(par.mask)])
            parm.est.region[, length(par.mask)] <- 
                exp(parm.est.region[, length(par.mask)])
        }
            
        z.pred.rnd <- matrix(NA, nrow(parm.est.rnd), sum(mask.test.rnd))
        z.pred.region <- matrix(NA, nrow(parm.est.region), 
                                sum(mask.test.region))
        for(k in 1 : nrow(parm.est.rnd)){
            z.pred.rnd[k, ] <- resp_pred(parm.est.rnd[k, ], 
                                         z.all[mask.train.rnd], 
                                         locxyz.all[mask.train.rnd, ], 
                                         locxyz.all[mask.test.rnd, ], 
                                         m, nuggets)$mu.pred
            z.pred.region[k, ] <- resp_pred(parm.est.region[k, ], 
                                            z.all[mask.train.region], 
                                            locxyz.all[mask.train.region, ], 
                                            locxyz.all[mask.test.region, ], 
                                            m, nuggets)$mu.pred
        }
        
        tbl[i, ] <- c(score_func(z.all[mask.test.rnd], z.pred.rnd), 
                      score_func(z.all[mask.test.region], z.pred.region))
        z.all.pred <- z.all
        z.all.pred[mask.test.rnd] <- colMeans(z.pred.rnd)
        if(plot.flag){
            sphere_plot(z.all.pred, rep(T, ns), grd.obj$lon, grd.obj$lat,
                        draw.palette = F, 
                        fn = paste0("rand-", type, "-", 
                                    mdl.lst[[i]], "-pred.pdf"))
        }else{
            input <- list(z.all = z.all.pred, mask = rep(T, ns), 
                          lon = grd.obj$lon, 
                          lat = grd.obj$lat, draw.palette = F, 
                          fn = paste0("rand-", type, "-", 
                                      mdl.lst[[i]], "-pred.pdf"))
            save(input, file = paste0("rand-", type, "-", 
                                      mdl.lst[[i]], "-pred.RData"))
        }
        
        z.all.pred <- z.all
        z.all.pred[mask.test.region] <- colMeans(z.pred.region)
        if(plot.flag){
            sphere_plot(z.all.pred, rep(T, ns), grd.obj$lon, grd.obj$lat,
                        draw.palette = F, 
                        fn = paste0("rect-", type, "-", mdl.lst[[i]], 
                                    "-pred.pdf"))
        }else{
            input <- list(z.all = z.all.pred, mask = rep(T, ns), 
                          lon = grd.obj$lon,
                          lat = grd.obj$lat, draw.palette = F, 
                          fn = paste0("rect-", type, "-", 
                                      mdl.lst[[i]], "-pred.pdf"))
            save(input, file = paste0("rect-", type, "-", 
                                      mdl.lst[[i]], "-pred.RData"))
        }
    }
    write.table(format(tbl, digits=3), paste0("tbl_", type, "_", seed, ".out"), 
                sep = " & ", eol = "\\\\\n", quote = F)
}
