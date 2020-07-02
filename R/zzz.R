#' @importFrom cli cat_line
#' @importFrom graphics abline axis box image layout lines par plot
#' @importFrom graphics points text title
#' @importFrom graphicsutils plot0 gpuPalette
#' @importFrom grDevices colorRampPalette dev.off png
#' @importFrom inSilecoMisc msgInfo msgSuccess scaleWithin keepLetters
#' @importFrom inSilecoMisc seqRg seqCol
#' @importFrom ks Hpi kde
#' @importFrom latex2exp TeX
#' @importFrom MASS lda
#' @importFrom mvtnorm dmvnorm rmvnorm
#' @importFrom orthomap orthomap
#' @importFrom progress progress_bar
#' @importFrom stats approx as.formula dchisq density dnorm integrate prcomp
#' @importFrom stats predict rnorm sd
#' @importFrom sf st_as_sf st_geometry
#' @importFrom utils combn

NULL

# HELPERS

output_dir <- function(dir = "output") {
  if (!dir.exists(dir)) {
    dir.create(dir)
    msgInfo("Folder", dir, "created!")
  }
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