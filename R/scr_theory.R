#' Theory script
#'
#' Code to reproduce analysis in the Supplementary information. Also, this
#' script creates figures S1, S2, S3 and S4.
#'
#' @param nrep number of replicates for figure S1, S3 and S4.
#' @param nrep2 number of replicates for figures S1, S5 and S6.
#' @param mx_sp_sz maximum sample size value.
#' @param rerun a logical. Should simulations be rerun for figures S5 and S6.
#'
#' @export

scr_theory <- function(nrep = 1e4, nrep2 = 1e5, mx_sp_sz = 50, rerun = FALSE) {

  # sequence of sample sizes
  nsz <- seq_len(mx_sp_sz)

  # create output dir if needed
  output_dir()

  ## Dissimilarity - MU
  msgInfo("Running Simulations for figure S1")
  res_sim <- list(
    twonorm_simu(nsz, nrep, c(0, .1), c(1, 1)),
    twonorm_simu(nsz, nrep, c(0, .25), c(1, 1)),
    twonorm_simu(nsz, nrep, c(0, .5), c(1, 1)),
    twonorm_simu(nsz, nrep, c(0, 1), c(1, 1)),
    twonorm_simu(nsz, nrep, c(0, 1), c(.5, .5))
  )
  res_ana <- list(
    mu_val(c(0, .1), 1, nsz),
    mu_val(c(0, .25), 1, nsz),
    mu_val(c(0, .5), 1, nsz),
    mu_val(c(0, 1), 1, nsz),
    mu_val(c(0, 1), .5, nsz)
  )
  plot_si("output/figS1.png", nsz, res_sim, res_ana)
  msgSuccess_fig("S1")


  ## Dissimilarity - SIGMA
  msgInfo("Creating figure S2")
  plot_log_ratio("output/figS2.png")
  msgSuccess_fig("S2")

  msgInfo("Running simulations for figure S3")
  res_sim <- list(
    twonorm_simu(nsz, nrep, c(0, 0), c(1, 1)),
    twonorm_simu(nsz, nrep, c(0, 0), c(1, 1.2)),
    twonorm_simu(nsz, nrep, c(0, 0), c(1, 1.5)),
    twonorm_simu(nsz, nrep, c(0, 0), c(1, 2)),
    twonorm_simu(nsz, nrep, c(0, 0), c(1, 5))
  )
  res_ana <- list(
    si_val(nsz, c(1, 1)),
    si_val(nsz, c(1, 1.2)),
    si_val(nsz, c(1, 1.5)),
    si_val(nsz, c(1, 2)),
    si_val(nsz, c(1, 5))
  )
  plot_si("output/figS3.png", nsz, res_sim, res_ana)
  msgSuccess_fig("S3")


  ## Dissimilarity - MU with 2, 5 and 10 observations
  msgInfo("Running simulations for figure S4")
  res_sim <- list(
    n_norm_simu(nsz, nrep, rep(0, 2), rep(1, 2), rep(.2, 2), rep(1, 2)),
    n_norm_simu(nsz, nrep, rep(0, 3), rep(1, 3), rep(.2, 3), rep(1, 3)),
    n_norm_simu(nsz, nrep, rep(0, 5), rep(1, 5), rep(.2, 5), rep(1, 5)),
    n_norm_simu(nsz, nrep, rep(0, 10), rep(1, 10), rep(.2, 10), rep(1, 10)),
    n_norm_simu(nsz, nrep, rep(0, 20), rep(1, 20), rep(.2, 20), rep(1, 20))
  )
  res_ana <- list(
    mu_val_n(c(0, 0), c(.2, .2), 1, nsz),
    mu_val_n(rep(0, 3), rep(.2, 3), 1, nsz),
    mu_val_n(rep(0, 5), rep(.2, 5), 1, nsz),
    mu_val_n(rep(0, 10), rep(.2, 10), 1, nsz),
    mu_val_n(rep(0, 20), rep(.2, 20), 1, nsz)
  )
  plot_si("output/figS4.png", nsz, res_sim, res_ana)
  msgSuccess_fig("S4")


  ## Correlation 2-10
  msgInfo("Running simulations for figure S5")
  scr_corr2("output/figS5.png", nrep2, rerun)
  msgSuccess_fig("S5")

  ## Correlation 1-10
  msgInfo("Running simulations for figure S6")
  if (rerun) {
    # takes ~3h for 1e5 repetitions (1 CPU i7)
    res <- simu_corr_10(nrep2)
    saveRDS(res, "inst/extdata/res_corr_10.rds")
  } else {
    res <- readRDS(system.file("extdata", "res_corr_10.rds",
      package = "spatialfingerprints"))
  }

  png("output/figS6.png", width = 5.5, height = 5, units = "in", res = 600)
    par(mar = c(4.5, 4.5, 1, 1), las = 1)
    plot(res$dim-res$cor, res$res, pch = 19, cex = .7, lwd = 2, lend = 1,
      ann = FALSE, axes = FALSE, ylim = c(.6, 1), yaxs = "i")
    axis(1)
    axis(2)
    box(bty = "l", lwd = 1.2)
    title(xlab = "Dimensionality", ylab = TeX("E(\\[A_1|S\\])"))
  dev.off()
  msgSuccess_fig("S6")

  invisible(NULL)
}




plot_si <- function(filename, seqx, res_sim, res_ana, cex_pt = 1.1) {
  png(filename, width = 5.5, height = 5, units = "in", res = 600)
    par(las = 1, bty = "l", mar = c(4.25, 4.25, 2, .5), mgp = c(2.6, .65, 0))
    n <- length(res_sim)
    pal <- colorRampPalette(c("black", "grey80"))(n)
    plot(range(seqx), y = c(.5, 1), cex = cex_pt, pch = 19, type = "n",
        xlab = "Sample size", ylab = TeX("$E(\\[A_1|S\\])$"))
    for (i in seq_len(n)) {
      points(seqx, res_sim[[i]], cex = cex_pt, pch = 19, col = pal[i])
      lines(seqx, res_ana[[i]], col = 2, lwd = 1.4)
    }
  dev.off()
}


plot_log_ratio <- function(filename) {
  n <- c(1, 2.5, 5, 10, 25)
  pal <- colorRampPalette(c("black", "grey80"))(length(n))
  sqx <- seq(0.01, 1, .01)
  sqxx <- c(rev(-sqx), 0, sqx)
  png(filename, width = 5.5, height = 5, units = "in", res = 600)
    par(las = 1, bty = "l", mar = c(4.5, 4.6, 2, .5), mgp = c(3.15, .65, 0),
      yaxs = "i", xaxs = "i", las = 1)

    plot(range(sqxx), c(.5, 1), type = "n", axes = FALSE, ylim = c(.49, 1.01),
      xlim = c(-1.04, 1.04), xlab = expression(frac(sigma[1], sigma[2])),
      ylab = TeX("$E(\\[A_1|S\\])$"))
    abline(v = 0, lty = 2)
    for (i in seq_along(n)) {
      lines(sqxx, log_ratios(sqx, n[i]), lwd = 1.4, col = pal[i])
    }
    axis(2)
    sql <- c(2, 5, 10)
    axis(1, at = c(rev(-log10(sql)), 0, log10(sql)),
      labels = c(rev(-sql), 1, sql))
    box(bty = "l")
  dev.off()
  invisible(NULL)
}


# Simulation fof figS5
scr_corr2 <- function(filename, nrep = 1e5, rerun = FALSE, cex_pt = .7) {

  if (rerun) {
    # takes ~1h for nrep = 1e5 (1CPU intel i7)
    df_res <- simu_corr_2(nrep)
    #  saveRDS(df_res, "inst/extdata/res_corr_2.rds")
  } else {
    df_res<- readRDS(
      system.file("extdata", "res_corr_2.rds", package = "spatialfingerprints")
    )
  }
  ls_res <- split(df_res$prob, f = df_res$corel)
  seqsize <- unique(df_res$sample_size)

  # create output dir if needed
  output_dir()

  png(filename, units = "in", width = 5.5, height = 5, res = 600)
    par(las = 1, mar = c(4.25, 4.25, 2, .5), mgp = c(2.6, .65, 0))
    pal2 <- rev(colorRampPalette(c("black", "grey70"))(length(ls_res)))
    plot0(c(0, 50), c(.5, 1))
    for (i in seq_along(ls_res)) {
      points(seqsize, ls_res[[i]], col = pal2[i], pch = 19, cex = cex_pt)
      lines(seqsize, ls_res[[i]], col = pal2[i], lwd = 1.4)
    }
    axis(1)
    axis(2)
    box(lwd = 1.2, bty = "l")
    title(xlab = "Sample Size", ylab = TeX("$E(\\[A_1|S\\])$"))
    text(10, 0.8, labels = TeX("$\\rho = 0.99$"), col = "black", pos = 4)
    text(10, 0.92, labels = TeX("$\\rho = 0$"), col = "grey70", pos = 2)
  dev.off()

  invisible(NULL)
}


simu_corr_10 <- function(nrep = 1e5, npt = 10, nvar = 10) {
  # correlation sequence
  seqc <- seq(0,.99, length = npt)
  #
  lsres <- list()
  ssz <- 25
  l <- 0
  d <- 1 # pic size
  Sigma <- matrix(.99, nvar, nvar)
  diag(Sigma) <- 1
  ##--
  mu1 <- rep(0, nvar)
  mu2 <- rep(1/sqrt(2*nvar), nvar)
  ##--
  pb <- progress_bar$new(format = paste0(cli::symbol$info,
      "computing [:bar] :percent eta: :eta"),
      total = (nvar - 1)*npt, clear = FALSE, width = 60)
  ##
  for (i in 2:nvar) {
    id <- (1:nvar)[-(2:i)]
    for (j in seq_along(seqc)) {
      pb$tick()
      Sigma[i,id] <- Sigma[id,i] <- 0.99 - seqc[j]
      tmp <- replicate(nrep, getProb(
          rmvnorm(ssz, mean = mu1, sigma = Sigma),
          mu1 = mu1,  mu2 = mu2, sigma = Sigma))
      l <- l+1
      ##-- build a data frame
      lsres[[l]] <- data.frame(
        dim = i,
        cor = 0.99 - seqc[j],
        res = mean(tmp),
        res_sd = sd(tmp)
        )
    }
  }
  do.call(rbind, lsres)
}




## Simulation helper functions

# vc_sz: vector of size
# nrep: number of replicates
# mu: mu parameters for the 2 distributions
# si: si parameters for the 2 distribution

twonorm_simu <- function(vc_sz, nrep, mu, si) {
  out <- list()
  for (i in seq_along(vc_sz)) {
    out[[i]] <- replicate(nrep, twonorm0(vc_sz[i], mu, si))
  }
  unlist(lapply(out, mean))
}

twonorm0 <- function(sz, mu, si) {
  r1 <- rnorm(sz, mu[1], si[1])
  1/(1 + exp(sum(
    dnorm(r1, mu[2], si[2], log = TRUE),
    - dnorm(r1, mu[1], si[1], log = TRUE)
    )))
}


# mu_1: vector of mu values (area A1 - true origin)
# si_1: vector of sigma values (area A1 - true origin)
# mu_2: vector of mu values (area A2)
# si_2: vector of sigma values (area A2)

n_norm_simu <- function(vc_sz, nrep, mu_1, mu_2, si_1, si_2) {
  nbio <- length(mu_1)
  stopifnot(all(lengths(list(si_1, mu_2, si_2)) == nbio))
  out <- list()
  for (i in seq_along(vc_sz)) {
    out[[i]] <- replicate(nrep, n_norm_0(vc_sz[i], mu_1, mu_2, si_1, si_2))
  }
  unlist(lapply(out, mean))
}


n_norm_0 <- function(sz, mu_1, si_1, mu_2, si_2) {

  mus_1 <- rep(mu_1, each = sz)
  sis_1 <- rep(si_1, each = sz)
  mus_2 <- rep(mu_2, each = sz)
  sis_2 <- rep(si_2, each = sz)
  r1 <- rnorm(length(mu_1)*sz, mus_1, sis_1)

  1/(1 + exp(sum(
    dnorm(r1, mus_2, sis_2, log = TRUE),
    - dnorm(r1, mus_1, sis_1, log = TRUE)
    )))
}


## Analytical results

# NB sigma (si) should be constant
mu_val <- function(mu, si, n) {
  unlist(lapply(n, function(x) logitnorm::momentsLogitnorm(
    x/(2 * si^2) * (mu[2]-mu[1])^2, sigma = sqrt(x) * abs(mu[2]-mu[1])/si)[1]))
}

# NB sigma (si) should be constant
mu_val_n <- function(mu_1, mu_2, si, n) {
  unlist(lapply(n, function(x) logitnorm::momentsLogitnorm(
    x/(2 * si^2) * sum((mu_2 - mu_1)^2), sigma = sqrt(x) * sqrt(sum(((mu_2 - mu_1)/si)^2)))[1]))
}

log_chi <- function(x, n, si) {
  ra <- si[1]/si[2]
  dchisq(x, df = n)*(1/(1 + (ra)^n*exp(- .5*((ra)*(ra) - 1)*x)))
}

si_val <- function(n, si) {
  unlist(
    lapply(n, function(x) integrate(log_chi, 0, Inf, n = x, si = si)$value)
  )
}


# function log(ratio), for figS3
log_ratios <- function(x, n) {
  x <- sort(x)
  vp <- vn <- double(length(x))
  for (i in seq_along(x)) {
    vp[i] <- integrate(log_chi, 0, Inf, n = n, si = c(1, 10^(x[i])))$value
    vn[i] <- integrate(log_chi, 0, Inf, n = n, si = c(1, 10^(-x[i])))$value
  }
  c(rev(vn), .5, vp)
}


# functions to be used with mutidimensional normal distribution

## LL for multivariate normal
LogLik2 <- function(val, mean, sigma) {
  sum(dmvnorm(val, mean, sigma = sigma, log = TRUE))
}

getProb <- function(val, mu1 = c(0,0), mu2 = c(1/sqrt(2), 1/sqrt(2)), sigma) {
    1/(1 + exp(LogLik2(val, mu2, sigma) - LogLik2(val, mu1, sigma)))
}




simu_corr_2 <- function(nrep = 1e4) {
  ### COVARIANCE
  nvar <- 2
  npt <- 6
  nss <- 10
  # correlation values
  seqc <- rev(seq(0, .99, length = npt))
  seqsize <- c(1:5, 10, 15, 20, 30, 40, 50)
  ##
  Sigma <- matrix(0, nvar, nvar)
  diag(Sigma) <- 1
  mu1 <- c(0, 0)
  ##
  df_res <- data.frame(
    sample_size = rep(seqsize, length(seqc)),
    corel = rep(seqc, each = length(seqsize)),
    prob = NA_real_,
    sd = NA_real_
  )

  pb <- progress_bar$new(format = paste0(cli::symbol$info,
      "computing [:bar] :percent eta: :eta"),
      total = nrow(df_res), clear = FALSE, width = 60)
  for (i in seq_len(nrow(df_res))) {
      pb$tick()
      ## correlation
      Sigma[cbind(c(1,2), c(2,1))] <- df_res$corel[i]
      ## use default mu values, i.e. c(0, 0) and c(1/sqrt(2), 1/sqrt(2))
      tmp <- replicate(nrep, getProb(
          rmvnorm(df_res$sample_size[i], mean = mu1, sigma = Sigma),
          sigma = Sigma))
      ##
      df_res$prob[i] <- mean(tmp)
      df_res$sd[i] <- sd(tmp)
  }

  df_res
}
