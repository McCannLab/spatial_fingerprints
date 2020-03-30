#' @importFrom graphics abline axis box lines par plot points
#' @importFrom graphicsutils plot0
#' @importFrom grDevices colorRampPalette dev.off png
#' @importFrom inSilecoMisc msgInfo msgSuccess
#' @importFrom latex2exp TeX
#' @importFrom mvtnorm dmvnorm rmvnorm
#' @importFrom progress progress_bar
#' @importFrom stats dchisq dnorm integrate rnorm

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