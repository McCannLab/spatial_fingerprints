#' Simulation functions
#'
#' @param nrep number of repetition for every iterations.
#' @param mxcb maximum number of biotracers combinations (see details).
#' @param nsample number of samples (individuals) to be tested (per region).
#' @param nbio number of biotracers to be used.
#' @param ndistr number of individuals used to generate distributions (per region).
#' @param noise noise to be added.
#'
#' @details
#' If the total number of combinations is smaller than `mxcb` then the
#' total number o of combinations is used.
#'
#' @export

scr_simu <- function() {
  # call all simulation
  invisible(0)
}

#' @describeIn scr_simu simulation for an increasing number of biotracers.
#' @export
simu_nbio <- function(nrep = 200, mxcb = 200, nsample, ndistr, noise) {
  out <- list(ksi = list(), lda = list())
  arg <- list(nsample = nsample, ndistr = ndistr, noise = noise)
  sq_nbio <- 1:17
  for (j in sq_nbio) {
    cat("--------------\nnbio = ", j, "\n")
    out$ksi[[j]] <- myreplic_combn(scr_ks_ind, arg, j,
        nrep = nrep, mxcb = mxcb, mxbio = 17, ngeo = 3)
    out$lda[[j]] <- myreplic_combn(scr_lda, arg, j,
        nrep = nrep, mxcb = mxcb, mxbio = 17, ngeo = 3)
  }
  out$nbio <- sq_nbio
  saveRDS(out, file = name_file("res_bio", nsample, ndistr, noise, nrep, mxcb))
  out
}

# simu_nbio(100, 100, 1, 27, 0.1)
# simu_nbio(500, 500, 1, 20, 0.5)
# simu_nbio(500, 500, 5, 20, 0)
# simu_nbio(500, 500, 5, 20, 0.5)


#' @describeIn scr_simu increase the number of samples used to build the distribution.
#' @export
simu_ndistr <- function(nrep = 200, mxcb = 200, nbio = 5, nsample, noise) {
  out <- list(ksi = list(), lda = list())
  sq_distr <- c(5:25)
  for (j in sq_distr) {
    cat("--------------\ndistr = ", j, "\n")
    arg <- list(nsample = nsample, ndistr = j, noise = noise)
    out$ksi[[j]] <- myreplic_combn(scr_ks_ind, arg, nbio,
        nrep = nrep, mxcb = mxcb, mxbio = 17, ngeo = 3)
    out$lda[[j]] <- myreplic_combn(scr_lda, arg, nbio,
        nrep = nrep, mxcb = mxcb, mxbio = 17, ngeo = 3)
  }
  out$ndistr <- sq_distr
  saveRDS(out, file = name_file("res_distr", nsample, "5-25", noise,
    nrep, mxcb))
  out
}
# simu_ndistr(100, 100, 10, 1, 0)
# simu_ndistr(500, 700, 3, 1, .5)
# simu_ndistr(500, 700, 3, 5, 0)
# simu_ndistr(500, 700, 3, 5, 0.5)
# simu_ndistr(5, 5, 5, 1, 0)
# res_dis <- readRDS("resdistr__ns1_ndvar_no0_nr100_mx100.rds") %>% get_res
# plot(c(5,25), c(1.5, 3), type = "n", ylab = "overall perf", xlab = "# distr")
# points( c(5, 10, 15, 20, 25), get_diagsum(res_dis$ksi), pch = 19, col = "gray10")
# points( c(5, 10, 15, 20, 25), get_diagsum(res_dis$lda), pch = 19, col = "gray60")


# nbio = 3 => 680 combination
# increase #sample used to infer

#' @describeIn scr_simu increase the number of samples used to build the distribution.
#' @export
simu_nsample <- function(nrep = 500, mxcb = 700, nbio = 3, ndistr, noise) {
  out <- list(ksi = list(), lda = list())
  sq_sample <- 1:10
  for (j in sq_sample) {
    cat("--------------\nnsample = ", j, "\n")
    arg <- list(nsample = j, ndistr = ndistr, noise = noise)
    out$ksi[[j]] <- myreplic_combn(scr_ks_ind, arg, nbio,
        nrep = nrep, mxcb = mxcb, mxbio = 17, ngeo = 3)
    out$lda[[j]] <- myreplic_combn(scr_lda, arg, nbio,
        nrep = nrep, mxcb = mxcb, mxbio = 17, ngeo = 3)
  }
  out$nsample <- sq_sample
  saveRDS(out, file = name_file("res_sample", "1-10", ndistr, noise,
    nrep, mxcb))
  out
}
# simul_nsample(500, 700, 3, 20, 0)
# simul_nsample(500, 700, 3, 20, .5)
# simul_nsample(500, 700, 3, 15, 0)
# simul_nsample(500, 700, 3, 15, 0.5)
# simul_nsample(100, 50, 5, 20, 0)
# res_dis <- get_res(readRDS("resdistr__nsvar_nd20_no0_nr100_mx50.rds"))
# plot(c(1,10), c(2, 3), type = "n", ylab = "overall perf", xlab = "# samples")
# points(1:10, get_diagsum(res_dis$ksi), pch = 19, col = "gray10")
# points(1:10, get_diagsum(res_dis$lda), pch = 19, col = "gray60")


# vc <- c(1:4, 6, 8, 10, 13, 16)
# out0 <- readRDS(file = "output/out_noise0_dis20_500_500.rds")
# out5 <- readRDS(file = "output/out_noise5_dis20_500_500.rds")
#
# res0 <- extract_res(out0)
# res5 <- extract_res(out5)
##
# par(mfrow = c(1, 2))
#
# pal <- c("gray10", "gray60")
# plot(c(0, 16), c(1.8, 3), xlab = "# Biotracers", ylab = "Overall performance (max = 3)", type = "n", main = "LDA")
# abline(h = 3, lty = 2, col = pal[2])
# points(vc, lapply(res0$lda, function(x) sum(diag(x))),
#   , pch = 19, col = pal[1])
# points(vc, lapply(res5$lda, function(x) sum(diag(x))), pc = 19, col = pal[2])
# legend("bottomright", c("Noise = 0", "Noise = .5", "#Sample = 1"), col = c(pal, NA), pch = 19, bty = "n", cex = 2)
#
#
# plot(c(0, 16), c(1.8, 3), xlab = "# Biotracers", ylab = "Overall performance (max = 3)", type = "n",main = "KSI")
# abline(h = 3, lty = 2, col = pal[2])
# points(vc, lapply(res0$ksi, function(x) sum(diag(x))),
#   , pch = 19, col = pal[1])
# points(vc, lapply(res5$ksi, function(x) sum(diag(x))), pc = 19, col = pal[2])


myreplic_combn <- function(FUN, arg, nbio, mxbio = 17, nrep = 1000, ngeo = 3,
  mxcb = 200) {
  tmp <- combn(mxbio, nbio)
  # seqCol(tmp)
  if (ncol(tmp) > mxcb) {
    tmp <- tmp[, sample(seqCol(tmp), mxcb), drop = FALSE]
  }
  out <- array(NA, c(ngeo, ngeo, ncol(tmp)))
  for (i in seqCol(tmp)) {
    print(i)
    arg$col_ids <- 2+tmp[,i]
    out[, ,i] <- myreplic(FUN, arg, nrep)
  }
  out
}

## where replicated
myreplic <- function(FUN, arg, nrep = 100) {
  apply(replicate(nrep, do.call(FUN, arg)), c(1,2), mean)
}

get_res <- function(x) {
  res <- list()
  res$ksi <- lapply(Filter(Negate(is.null), x$ksi), apply, c(1,2), mean)
  res$lda <- lapply(Filter(Negate(is.null), x$lda), apply, c(1,2), mean)
  res
}

get_diagsum <- function(x) unlist(lapply(x, function(x) sum(diag(x))))
