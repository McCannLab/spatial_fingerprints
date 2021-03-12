#' Analysis pipeline
#'
#' Run all analyses. Note this requires the main (and time-consuming) analyses to be run first. 
#'
#' @export

pipeline <- function() {
  #
  set.seed(7891)
  msgInfo("Format data...")
  
  # centered and scaled (cs)
  write.csv(get_data_ready(), file = "inst/extdata/data_f_cs.csv",
  row.names = FALSE)
  
  # apply a pca
  write.csv(get_data_ready(pca = TRUE), file = "inst/extdata/data_f_pca.csv",
  row.names = FALSE)
  # format data
  
  #
  df_all <- system.file("extdata", "isoready.csv", package =
    "spatialfingerprints")

  msgInfo("Main Figures...")
  scr_fig_concept()
  scr_fig2()
  scr_fig3()
  scr_fig4()
  scr_fig5()


  msgInfo("SI Figures...")
  scr_theory()
  scr_figS6()
  #
  invisible(NULL)
}