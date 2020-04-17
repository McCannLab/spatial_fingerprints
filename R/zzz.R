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
#' @importFrom stats as.formula dchisq density dnorm integrate predict rnorm sd
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
