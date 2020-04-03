#' Analysis pipeline
#'
#' Run all analyses.
#'
#' @export

pipeline <- function() {
  #
  set.seed(7891)
  df_all <- system.file("extdata", "isoready.csv", package =
    "spatialfingerprints")
  #
  msgInfo("Main Figures...")
  script_fig_concept()


  msgInfo("SI Figures...")
  scr_theory()
  #
  invisible(NULL)
}