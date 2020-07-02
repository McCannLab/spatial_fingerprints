#' Simulation functions
#'
#' @param method statistical approach, either "lda", "nb" or "ml".
#' @param nrep number of repetition for every iterations.
#' @param mxcb maximum number of biotracers combinations (see details).
#' @param nsample number of samples (individuals) to be tested (per region).
#' @param ndistr number of individuals used to generate distributions (per region).
#' @param noise noise to be added.
#'
#' @details
#' If the total number of combinations is smaller than `mxcb` then the
#' total number of combinations is used.
#'
#' @export
#' @example
#' simu_nbio("lda")


simu_nbio <- function(method = c("lda", "nb", "ml"), nrep = 100, mxcb = 10,
  nsample = 10, ndistr = 20, noise = 0) {

  df_dat <- get_data_ready()
  method <- match.arg(method)
  arg <- list(method = method, nsample = nsample, ndistr = ndistr,
      noise = noise, df_dat = df_dat)
  sq_nbio <- seq_len(17)
  out <- list()
  for (j in sq_nbio) {
      cat_line(cli::symbol$star, " nbio = ", j)
      # average over nrep replicates
      out[[j]] <- myreplic_combn(find_origin, arg, j, nrep = nrep, mxcb = mxcb,
          mxbio = 17, ngeo = 3)
  }
  out
}






# helpers

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
