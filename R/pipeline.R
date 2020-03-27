#' Analysis pipeline
#'
#' Run all analyses.
#'
#' @export

pipeline <- function() {
  #
  set.seed(7891)
  #
  msgInfo("SI Figures")
  scr_theory()
  #
  invisible(NULL)
}