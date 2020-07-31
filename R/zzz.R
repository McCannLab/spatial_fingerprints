#' @importFrom cli cat_line style_underline
#' @importFrom graphics abline axis box image layout lines par plot
#' @importFrom graphics points text title legend boxplot mtext
#' @importFrom graphicsutils plot0 gpuPalette darken
#' @importFrom grDevices colorRampPalette dev.off png
#' @importFrom inSilecoMisc msgInfo msgSuccess scaleWithin keepLetters
#' @importFrom inSilecoMisc seqRg seqCol adjustStrings
#' @importFrom ks Hpi kde
#' @importFrom latex2exp TeX
#' @importFrom MASS lda
#' @importFrom mvtnorm dmvnorm rmvnorm
#' @importFrom orthomap orthomap
#' @importFrom progress progress_bar
#' @importFrom stats approx as.formula dchisq density dnorm integrate prcomp
#' @importFrom stats aggregate predict rnorm sd coefficients deviance
#' @importFrom stats dist lm nls
#' @importFrom sf st_as_sf st_geometry
#' @importFrom utils combn write.csv


NULL



# HELPERS

output_dir <- function(dir = "output", recursive = TRUE) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = recursive)
    msgInfo("Folder", style_underline(dir), "created!")
  } else msgInfo("Folder", style_underline(dir), "has already been created!")
  invisible(dir)
}

msgSuccess_fig <- function(n, dir = "output")
  msgSuccess("Fig", n, "created!", paste0("See ", dir, "/fig", n, ".png"))


name_file <- function(pref = "res", nsample, ndistr, noise, nrep, mxcb) {
  paste0(pref, "_nsa", nsample, "_ndi", ndistr, "_noi", noise,
    "_nre", nrep, "_mxc", mxcb, ".rds")
}


# sampling combination of biotracer helper
get_replic <- function(n, p, mx = 200) {
  tmp <- choose(n, p)
  if (mx <= tmp) {
    combn(n, p)[, sample(tmp, mx), drop = FALSE]
  } else combn(n, p)
}



##
# add noise
add_noise <- function(train, mean = 0, sd = 1) {
  apply(train, 2, function(x) x + rnorm(length(x), mean, sd))
}

add_noise_v <- function(train, mean = 0, sd = 1) {
  train + rnorm(length(train), mean, sd)
}


# convert a matrix of -log likelihood into a proba for each state
# NB: FOR EACH ELEMENT row i column j
# $exp(-m_{i,j})/\sum_j{exp(-m_{i,j})}$
toprob <- function(x) {
  tmp <- exp(-x)
  # NB: call to t() required, otherwise results are in the wrong order
  t(apply(tmp, 1, function(x) x/sum(x)))
}


julia_call <- function(pca, method, n_bio, n_distr, n_sample, mxcb, n_rep, inp, out = "output/res_ml/") {
  output_dir(out)
  inp <- ifelse(pca,
      "inst/exdata/data_f_pca.csv",
      "inst/exdata/data_f_cs.csv"
    )
  call <- paste('julia -p 1 scr/main.jl', pca * 1, n_bio, n_distr, n_sample, mxcb, n_rep, inp, out)
  msgInfo('Now computing via julia scripts')
  system(call)
  invisible(call)
}







## Helper functions scr_fig5()

seq_rg <- function(x, n, offset = 0) {
  rg <- range(x)
  seq(rg[1] - offset, rg[2] + offset, length.out = n)
}

myvar <- function(x) sum((x - mean(x))^2)


## max perf
getmax <- function(..., perf, off = -2) {
  ids <- list(...)
  max(perf[unlist(ids) + off])
}
## sum perf
getsum <- function(..., perf, off = -2) {
  ids <- list(...)
  sum(perf[unlist(ids) + off])
}
# isotops / FA or both?
get_categ2 <- function(..., off = -2) {
  ids <- list(...)
  sum(unlist(ids) %in% 1:3)
}

inertia_tot <- function(x, y) {
  sum((x - matrix(1, nrow = nrow(x)) %*% y)^2)
}

mean_dist_cent_wtn <- function(x, y) {
  apply((x - matrix(1, nrow = nrow(x)) %*% y)^2, 2, mean)
}

inertia <- function(x, grp = rep(1:3, each = 30)) {
  cent_all <- apply(x, 2, mean)
  inert_tot <- inertia_tot(x, cent_all)

  x_grp <- split(x, grp)
  ng <- length(x_grp)
  cent_grp <- lapply(x_grp, apply, 2, mean)

  inert_wtn_grp <- mapply(FUN = inertia_tot, x_grp, cent_grp)
  inert_wtn <- sum(inert_wtn_grp)

  # I do the following so as to double check wtn/btw otherwise I could have used tot-wtn
  inert_btw_grp <- unlist(lapply(x_grp, nrow)) *
    unlist(lapply(cent_grp, function(x) sum((x - cent_all)^2)))
  inert_btw <- sum(inert_btw_grp)


  # mean dist between centroids
  mean_dist_cent <- t(mapply(mean_dist_cent_wtn, x_grp, cent_grp))
  dist_cent <- mean(dist(do.call(rbind, cent_grp))^2)

  # compute mean d^2 in all dimension for all pair of centroid minus the mean d^2 of all distance for all dimension for the 2 regions
  res <- list()
  k = 0
  for (i in 1:(ng-1)) {
    for (j in (i+1):ng) {
      k = k+1
      res[[k]] <- (cent_grp[[i]] - cent_grp[[j]])^2 - mean_dist_cent[i,] - mean_dist_cent[j,]
    }
  }
  dist_cent_cor <- sum(apply(do.call(rbind, res), 2, mean))
  c(
    inert_tot = inert_tot ,
    inert_wtn = inert_wtn,
    inert_btw = inert_btw,
    dist_cent = dist_cent,
    dist_cent_cor = dist_cent_cor
  )
}

getdist <- function(hh) {
  m1 <- apply(hh[1:30, ], 2, mean)
  m2 <- apply(hh[31:60, ], 2, mean)
  m3 <- apply(hh[61:90, ], 2, mean)
  mean(sum((m1-m2)^2), sum((m1-m3)^2), sum((m2-m3)^2))
}

fitexp <- function(data, ...) {
  mod <- nls(perf ~ 1 - .66*exp(-d*mean_distb), list(d = 1), data = data)
  d <- coefficients(mod)
  x <- seq_rg(data$mean_distb, 100)
  lines(x, 1 - .66*exp(-d*x), ...)
  invisible(mod)
}

fitexp2 <- function(data, cl, ...) {
  mod <- nls(perf ~ 1-.66 * exp(- d*cl), list(d = 1), data = data)
  d <- coefficients(mod)
  x <- seq_rg(cl, 100)
  lines(x, 1 - .66*exp(-d*x), ...)
  invisible(mod)
}