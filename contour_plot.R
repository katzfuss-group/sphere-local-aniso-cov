rm(list = ls())
devtools::install_github("https://github.com/katzfuss-group/GPvecchia",
  ref = "JJ"
)
source("func.R")

plot_contour <- function(par, locs, level = 0.7, fn = NULL) {
  n_locs <- nrow(locs)
  contours <- list()
  if (!is.null(fn)) {
    fig.width <- 13.44
    fig.height <- 8.73
    pdf(fn, fig.width, fig.height)
  }

  mapPlot(
    longitude = c(-180, 180, 0, 0), latitude = c(0, 0, -90, 90),
    grid = TRUE, col = "lightgray", drawBox = FALSE,
    longitudelim = c(-180, 180), type = "n"
  )
  for (i in 1:n_locs) {
    lon <- locs[i, 1]
    lat <- locs[i, 2]
    lon_range <- seq(from = lon - 0.5, to = lon + 0.5, length.out = 100)
    lat_range <- seq(from = lat - 0.5, to = lat + 0.5, length.out = 100)
    lon_range <- lon_range[lon_range < pi & lon_range > -pi]
    lat_range <- lat_range[lat_range < pi / 2 & lat_range > -pi / 2]
    grid <- as.matrix(expand.grid(lon_range, lat_range))
    grid <- rbind(locs[i, , drop = F], grid)
    alpha <- par[1:3]
    beta <- par[4:6]
    kappa <- par[7]
    nu <- par[8]
    range <- par[9]
    cov_mat <- cov_gen_row1(grid, alpha, beta, kappa, nu, range, 0)
    idx <- cov_mat[1, ] < level + 0.01 & cov_mat[1, ] > level - 0.01
    if (sum(idx) == 0) {
      lon_range <- seq(from = lon - 0.1, to = lon + 0.1, length.out = 100)
      lat_range <- seq(from = lat - 0.1, to = lat + 0.1, length.out = 100)
      lon_range <- lon_range[lon_range < pi & lon_range > -pi]
      lat_range <- lat_range[lat_range < pi / 2 & lat_range > -pi / 2]
      grid <- as.matrix(expand.grid(lon_range, lat_range))
      grid <- rbind(locs[i, , drop = F], grid)
      cov_mat <- cov_gen_row1(grid, alpha, beta, kappa, nu, range, 0)
      idx <- cov_mat[1, ] < level + 0.02 & cov_mat[1, ] > level - 0.02
    }
    if (sum(idx) <= 1) {
      stop(cat("Cannot find contour of level", level))
    }
    grid_sub <- grid[idx, , drop = F]
    lon_sub <- sort(unique(grid_sub[, 1]))
    grid_sub_tmp <- matrix(NA, 0, 2)
    for (lon_tmp in lon_sub) {
      idx <- grid_sub[, 1] == lon_tmp & grid_sub[, 2] >= lat
      if (sum(idx) > 0) {
        grid_sub_tmp <- rbind(grid_sub_tmp, c(lon_tmp, max(grid_sub[idx, 2])))
      }
    }
    for (lon_tmp in rev(lon_sub)) {
      idx <- grid_sub[, 1] == lon_tmp & grid_sub[, 2] < lat
      if (sum(idx) > 0) {
        grid_sub_tmp <- rbind(grid_sub_tmp, c(lon_tmp, min(grid_sub[idx, 2])))
      }
    }
    grid_sub_tmp <- rbind(grid_sub_tmp, grid_sub_tmp[1, ])
    # mapPoints(lon * 180 / pi, lat * 180 / pi, col = "red", pch = 16)
    mapLines(grid_sub_tmp[, 1] * 180 / pi, grid_sub_tmp[, 2] * 180 / pi, col = "red")
  }
  if (!is.null(fn)) {
    dev.off()
  }
}

cov_gen_row1 <- function(grd.all, alpha, beta, kappa, nu, range, nuggets,
                         cluster = NULL) {
  n <- nrow(grd.all)
  gamma1 <- exp(alpha[1] + alpha[2] * sin(grd.all[, 1]) + alpha[3] * grd.all[, 2])
  gamma2 <- exp(beta[1] + beta[2] * sin(grd.all[, 1]) + beta[3] * grd.all[, 2])
  if (is.null(cluster)) {
    inner.cluster <- T
    cluster <- makeCluster(detectCores())
  }

  Sigma.input <- cbind(grd.all, gamma1, gamma2, kappa)
  Sigma.input <- lapply(
    1:nrow(Sigma.input),
    function(x) {
      return(Sigma.input[x, ])
    }
  )
  Sigma.func <- function(x) {
    D <- diag(c(1, x[3], x[4]))
    Sigmatilde <- Rx(x[5]) %*% D %*% t(Rx(x[5]))
    Rz(x[1]) %*% Ry(x[2]) %*% Sigmatilde %*% t(Ry(x[2])) %*%
      t(Rz(x[1]))
  }
  clusterExport(cluster, c("Rx", "Ry", "Rz", "cart"), envir = .GlobalEnv)
  Sigma <- lapply(Sigma.input, Sigma.func)
  clusterExport(cluster, c("grd.all", "Sigma", "n"), envir = environment())
  dis.func <- function(x) {
    if (x < 1) {
      return(0)
    }
    i <- (x - 1) %% n + 1
    j <- floor((x - 1) / n) + 1
    cart.diff <- cart(grd.all[i, 1], grd.all[i, 2]) -
      cart(grd.all[j, 1], grd.all[j, 2])
    sqrt(2) * (cart.diff %*% solve(Sigma[[i]] + Sigma[[j]]) %*%
      t(cart.diff))^(1 / 2)
  }
  c.func <- function(x) {
    if (x < 1) {
      return(0)
    }
    i <- (x - 1) %% n + 1
    j <- floor((x - 1) / n) + 1
    (det(Sigma[[i]]))^(1 / 4) * (det(Sigma[[j]]))^(1 / 4) *
      (det((Sigma[[i]] + Sigma[[j]]) / 2))^(-1 / 2)
  }
  tmp.mat <- matrix(1:n, n, n, F) + outer(rep(n, n), 0:(n - 1))
  tmp.mat[lower.tri(tmp.mat)] <- 0
  tmp.mat <- tmp.mat[1, , drop = F]
  dis.mat <- apply(
    X = tmp.mat, MARGIN = c(1, 2),
    FUN = dis.func
  )
  c.mat <- apply(
    X = tmp.mat, MARGIN = c(1, 2),
    FUN = c.func
  )
  C <- c.mat * Matern(dis.mat, range = range, nu = nu)
  C
}
locs <- as.matrix(expand.grid(
  c(-2:2) * pi / 3,
  c(-3:3) * pi / 8
))

par_iso <- c(8.476336, 0, 0, 8.476336, 0, 0, 0, 2.5, exp(0.07437005))
par_axial <- c(
  -2.628519, 0, -0.6661641, -6.35629, 0, -0.2518264, 0, 2.5,
  exp(0.07437005)
)
par_nonsta <- c(
  -2.441852, -0.0654426, -0.7878552, -6.375738, -0.2483788,
  -0.2730235, -6.32839, 2.5, exp(0.07437005)
)
for (type in c("iso", "axial", "nonsta"))
# for(type in c("iso"))
{
  plot_contour(get(paste0("par_", type)), locs, fn = paste0("contour_", type, ".pdf"))
}
