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
  nsz = seq_len(mx_sp_sz)

  ## mu
  msgInfo("Running Simulation for figure S1")
  res1 <- twonorm_simu(nsz, nrep, c(0, .5), c(2, 2))
  res2 <- twonorm_simu(nsz, nrep, c(0, .5), c(1, 1))
  res3 <- twonorm_simu(nsz, nrep, c(0, 1), c(1, 1))

  # create output dir if need
  output_dir()

  png("output/fig_S1.png", width = 7, height = 6, units = "in", res = 300)
    par(las = 1, bty = "l", mar = c(4.5, 4.5, 2, .5))
    plot(nsz, res1, ylim = c(.5, 1), cex = 2, pch = 19,
      xlab = "Sample size", ylab = TeX("$E(\\[A_1|S\\])$"))
    lines(nsz, mu_val(c(0, .5), 2 , nsz), col = 2, lwd = 2)
  #
    points(nsz, res2, cex = 2, pch = 19, col = "grey30")
    lines(nsz, mu_val(c(0, .5), 1 , nsz), col = 2, lwd = 2)
  #
    points(nsz, res3, cex = 2, pch = 19, col = "grey60")
    lines(nsz, mu_val(c(0, 1), 1 , nsz), col = 2, lwd = 2)
  dev.off()

  msgSuccess_fig("S1")


  ## si
  msgInfo("Running Simulation for figure S2")
  res1 <- twonorm_simu(nsz, nrep, c(0, 0), c(1, 1.2))
  res2 <- twonorm_simu(nsz, nrep, c(0, 0), c(1, 1.5))
  res3 <- twonorm_simu(nsz, nrep, c(0, 0), c(1, 2))

  png("output/fig_S2.png", width = 7, height = 6, units = "in", res = 300)
    par(las = 1, bty = "l", mar = c(4.5, 4.5, 2, .5))
    plot(nsz, res1, ylim = c(.5, 1), cex = 2, pch = 19,
        xlab = "Sample size", ylab = TeX("$E(\\[A_1|S\\])$"))
    lines(nsz, si_val(nsz, c(1, 1.2)), col = 2, lwd = 2)
  #
    points(nsz, res2, cex = 2, pch = 19, col = "grey30")
    lines(nsz, si_val(nsz, c(1, 1.5)), col = 2, lwd = 2)
  #
    points(nsz, res3, cex = 2, pch = 19, col = "grey60")
    lines(nsz, si_val(nsz, c(1, 2)), col = 2, lwd = 2)
  dev.off()

  msgSuccess_fig("S2")

  invisible(NULL)
}



mu_val <- function(mu, si, n) {
  unlist(lapply(n, function(x) logitnorm::momentsLogitnorm(
    x/(2*si^2)*(mu[2]-mu[1])^2, sigma = sqrt(x)*abs(mu[2]-mu[1])/si)[1]))
}

twonorm_simu <- function(val, nrep, mu, si) {
  out <- list()
  for (i in seq_along(val)) {
    out[[i]] <- replicate(nrep, twonorm0(val[i], mu, si))
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

log_chi <- function(x, n, si) {
  ra <- si[1]/si[2]
  dchisq(x, df = n)*(1/(1 + (ra)^n*exp(- .5*((ra)*(ra) - 1)*x)))
}

si_val <- function(n, si) {
  unlist(
    lapply(n, function(x) integrate(log_chi, 0, Inf, n = x, si = si)$value)
  )
}