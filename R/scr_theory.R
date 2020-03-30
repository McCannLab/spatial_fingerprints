#' Theory script
#'
#' Code to reproduce analysis in the Supplementary information. Also, this
#' script creates figures S1, S2, S3 and S4.
#'
#' @param nrep number of replicates.
#' @param mx_sp_sz maximum sample size value.
#'
#' @export

scr_theory <- function(nrep = 1e4, mx_sp_sz = 50) {

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
    n_norm_simu(nsz, nrep, rep(0, 2), rep(1, 2), c(.2, .5), rep(1, 2)),
    n_norm_simu(nsz, nrep, rep(0, 20), rep(1, 20), rep(.2, 20), rep(1, 20))
  )
  res_ana <- list(
    mu_val_n(c(0, 0), c(.2, .2), 1, nsz),
    mu_val_n(rep(0, 3), rep(.2, 3), 1, nsz),
    mu_val_n(rep(0, 5), rep(.2, 5), 1, nsz),
    mu_val_n(c(0, 0), c(.2, .5), 1, nsz),
    mu_val_n(rep(0, 20), rep(.2, 20), 1, nsz)
  )
  plot_si("output/figS4.png", nsz, res_sim, res_ana)
  msgSuccess_fig("S4")


  ## Correlation
  msgInfo("Running simulations for figure S5")



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





# integrate(log_chi, 0, Inf, n = x, si = si)$value)
#
#
# foo <- function(x, r, n) {
#   1/(1 + r^(n/2) * exp(-.5 * (r-1)*x))
# }


# sqx <- 10^seq(-2, 2, .01)
# plot(sqx, foo(sqx, 2, 10), type = "l")
# lines(sqx, foo(sqx, .2, 10), type = "l")
# lines(sqx, foo(sqx, 1, 10), type = "l")



# ```
# from sage.functions.other import symbolic_limit as slimit
# n = var('n')
# x = var('x')
# y = var('y')
# f(y) = (1/2)^(n/2)*(1/gamma(n/2, hold = True))*y^(n/2-1)*exp(-y/2)/(1+x^n*exp(-1/2*(x^2-1)*y))
# assume(n = 1)
# assume(x > 1)
# integrate(f(y), y, 0, infinity)
# ```