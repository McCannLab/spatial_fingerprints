#' Simulation functions
#'
#' @param method statistical approach, either "lda", "nb" or "ml".
#' @param nrep number of repetition for every iterations.
#' @param mxcb maximum number of biotracers combinations (see details).
#' @param nsample number of samples (individuals) to be tested (per region).
#' @param ndistr number of individuals used to generate distributions (per region).
#' @param nbio number of biotracers to be used.
#' @param noise noise to be added.
#' @param pca a logical. See [get_data_ready()].
#'
#' @details
#' If the total number of combinations is smaller than `mxcb` then the
#' total number of combinations is used.
#'
#' @export
#' @examples
#' \dontrun{
#' simu_nbio("lda")
#' simu_ndistr("lda")
#' simu_nsample()
#' }

simu_nbio <- function(method = c("lda", "nb", "ml"), nrep = 100, mxcb = 10,
  nsample = 10, ndistr = 20, noise = 0, pca = FALSE) {

  method <- match.arg(method)
  df_dat <- get_data_ready(pca = pca)
  out <- list()
  arg <- list(method = method, nsample = nsample, ndistr = ndistr,
      noise = noise, df_dat = df_dat)

  sq_nbio <- seq_len(17)
  for (j in sq_nbio) {
      cat_line(cli::symbol$star, " nbio = ", j)
      # average over nrep replicates
      out[[j]] <- myreplic_combn(find_origin, arg, j, nrep = nrep, mxcb = mxcb,
          mxbio = 17, ngeo = 3)
  }

  out
}

#' @describeIn simu_nbio same as [simu_nbio()] but use pca and keep the axes ordered.
#' @export
simu_nbio_order <- function(method = c("lda", "nb", "ml"), nrep = 100,
  nsample = 10, ndistr = 20, noise = 0) {

  method <- match.arg(method)
  df_dat <- get_data_ready(pca = TRUE)
  out <- list()
  arg <- list(method = method, nsample = nsample, ndistr = ndistr,
      noise = noise, df_dat = df_dat)

  sq_nbio <- seq_len(17)
  for (j in sq_nbio) {
      cat_line(cli::symbol$star, " nbio = ", j)
      arg <- list(method = method, df_dat = df_dat, nsample = nsample,
          ndistr = ndistr, col_ids = seq_len(j) + 2, noise = noise)
      out[[j]] <- replicate(nrep, do.call(find_origin, arg))
  }

  out
}



#' @describeIn simu_nbio effect of the number of samples used to build the distribution.
#' @export
simu_ndistr <- function(method = c("lda", "nb", "ml"), nrep = 20, mxcb = 20,
  nbio = 5, nsample = 1, noise = 0, pca = FALSE) {

  method <- match.arg(method)
  df_dat <- get_data_ready(pca = pca)
  out <- list()

  sq_distr <- 5:25
  for (j in sq_distr) {
    cat_line(cli::symbol$star, " ndistr = ", j)
    arg <- list(method = method, df_dat = df_dat, nsample = nsample, ndistr = j, noise = noise)
    out[[j]] <- myreplic_combn(find_origin, arg, nrep = nrep, nbio, mxcb = mxcb,
        mxbio = 17, ngeo = 3)
  }

  out
}

#' @describeIn simu_nbio effect of the number of samples used for inference.
simu_nsample <- function(method = c("lda", "nb", "ml"), nrep = 20, mxcb = 20, nbio = 5, ndistr = 20, noise = 0, pca = FALSE) {

  method <- match.arg(method)
  df_dat <- get_data_ready(pca = pca)
  out <- list()

  sq_sample <- 1:10
  for (j in sq_sample) {
    cat_line(cli::symbol$star, " nsample = ", j)
    arg <- list(method = method, df_dat = df_dat, nsample = j, ndistr = ndistr,
        noise = noise)
    out[[j]] <- myreplic_combn(find_origin, arg, nbio, nrep = nrep, mxcb = mxcb,
        mxbio = 17, ngeo = 3)
  }

  out
}

#' @describeIn simu_nbio effect of the number of samples used for inference.
simu_noise <- function(method = c("lda", "nb", "ml"), nrep = 20, mxcb = 20, nbio = 5, ndistr = 20, nsample = 5, pca = FALSE) {

  method <- match.arg(method)
  df_dat <- get_data_ready(pca = pca)
  out <- list()

  sq_noise <- 10^seq(-4, 0, .25)
  for (j in seq_along(sq_noise)) {
    cat_line(cli::symbol$star, " noise = ", sq_noise[j])
    arg <- list(method = method, df_dat = df_dat, nsample = nsample,
        ndistr = ndistr, noise = sq_noise[j])
    out[[j]] <- myreplic_combn(find_origin, arg, nbio, nrep = nrep, mxcb = mxcb,
        mxbio = 17, ngeo = 3)
  }

  out
}




# helpers

# replicate for a given combination
myreplic_combn <- function(FUN, arg, nbio, mxbio = 17, nrep = 1000, ngeo = 3,
  mxcb = 200) {
  tmp <- get_replic(mxbio, nbio, mxcb)
  out <- array(NA, c(ngeo, ngeo, ncol(tmp)))
  for (i in seqCol(tmp)) {
    arg$col_ids <- 2 + tmp[, i]
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
