#' Theory script
#'
#' Code to reproduce analysis for figures S1 and S2.
#'
#' @param nrep number of replicates.
#' @param mx_sp_sz maximum sample size value.
#'
#' @export
scr_theory <- function(nrep = 1e4, mx_sp_sz = 50) {

  # sequence of sample sizes
  nsz <- seq_len(mx_sp_sz)

  cex_pt <- 1.25

  ## mu
  msgInfo("Running Simulation for figure S1")
  res1 <- twonorm_simu(nsz, nrep, c(0, .5), c(2, 2))
  res2 <- twonorm_simu(nsz, nrep, c(0, .5), c(1, 1))
  res3 <- twonorm_simu(nsz, nrep, c(0, 1), c(1, 1))

  # create output dir if needed
  output_dir()

  ## Dissimilarity - MU
  png("output/fig_S1.png", width = 5.5, height = 5, units = "in", res = 300)
    par(las = 1, bty = "l", mar = c(4.5, 4.5, 2, .5))
    plot(nsz, res1, ylim = c(.5, 1), cex = cex_pt, pch = 19,
      xlab = "Sample size", ylab = TeX("$E(\\[A_1|S\\])$"))
    lines(nsz, mu_val(c(0, .5), 2 , nsz), col = 2, lwd = 2)
  #
    points(nsz, res2, cex = cex_pt, pch = 19, col = "grey30")
    lines(nsz, mu_val(c(0, .5), 1 , nsz), col = 2, lwd = 2)
  #
    points(nsz, res3, cex = cex_pt, pch = 19, col = "grey60")
    lines(nsz, mu_val(c(0, 1), 1 , nsz), col = 2, lwd = 2)
  dev.off()

  msgSuccess_fig("S1")


  ## Dissimilarity - SIGMA
  msgInfo("Running Simulation for figure S2")
  res1 <- twonorm_simu(nsz, nrep, c(0, 0), c(1, 1.2))
  res2 <- twonorm_simu(nsz, nrep, c(0, 0), c(1, 1.5))
  res3 <- twonorm_simu(nsz, nrep, c(0, 0), c(1, 2))

  png("output/fig_S2.png", width = 5.5, height = 5, units = "in", res = 300)
    par(las = 1, bty = "l", mar = c(4.5, 4.5, 2, .5))
    plot(nsz, res1, ylim = c(.5, 1), cex = cex_pt, pch = 19,
        xlab = "Sample size", ylab = TeX("$E(\\[A_1|S\\])$"))
    lines(nsz, si_val(nsz, c(1, 1.2)), col = 2, lwd = 2)
  #
    points(nsz, res2, cex = cex_pt, pch = 19, col = "grey30")
    lines(nsz, si_val(nsz, c(1, 1.5)), col = 2, lwd = 2)
  #
    points(nsz, res3, cex = cex_pt, pch = 19, col = "grey60")
    lines(nsz, si_val(nsz, c(1, 2)), col = 2, lwd = 2)
  dev.off()

  msgSuccess_fig("S2")


  ## Dissimilarity - MU with 2, 5 and 10 observations
  msgInfo("Running Simulation for figure S3")
  res1 <- n_norm_simu(nsz, nrep, rep(0, 2), rep(1, 2), rep(.2, 2), rep(1, 2))
  res2 <- n_norm_simu(nsz, nrep, rep(0, 2), rep(1, 2), c(.2, .4), rep(1, 2))
  res3 <- n_norm_simu(nsz, nrep, rep(0, 10), rep(1, 10), rep(.2, 10), rep(1, 10))

  png("output/fig_S3.png", width = 5.5, height = 5, units = "in", res = 600)
    par(las = 1, bty = "l", mar = c(4.5, 4.5, 2, .5))
    plot(nsz, res1, ylim = c(.5, 1), cex = cex_pt, pch = 19,
        xlab = "Sample size", ylab = TeX("$E(\\[A_1|S\\])$"))
    lines(nsz, mu_val_n(c(0, 0), c(.2, .2),  1, nsz), col = 2, lwd = 2)
  #
    points(nsz, res2, cex = cex_pt, pch = 19, col = "grey30")
    lines(nsz, mu_val_n(c(0, 0), c(.2, .4),  1, nsz), col = 2, lwd = 2)
    # lines(nsz, mu_val_n(rep(0, 5), rep(.2, 5),  1, nsz), col = 2, lwd = 2)
  #
    points(nsz, res3, cex = cex_pt, pch = 19, col = "grey60")
    lines(nsz, mu_val_n(rep(0, 10), rep(.2, 10),  1, nsz), col = 2, lwd = 2)
  dev.off()


  msgSuccess_fig("S3")
  invisible(NULL)
}











## Simulation functions

# vc_sz: vector of size
# nrep: number of replicates
# mu: mu paramters for the 2 distributions
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

