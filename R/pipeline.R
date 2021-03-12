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
  
  # step that cannot be reproduced due to the absence of data were 
  # commented out.
  # Note that in the code there are comments to explain how SI figures were 
  # reproduced.
  msgInfo("Main Figures...")
  scr_fig_concept()
  # scr_fig2()
  # scr_fig3()
  # scr_fig4()
  # scr_fig5()

  msgInfo("SI Figures...")
  scr_theory()
  # scr_figS6()
  # scr_figS7()
  #
  invisible(NULL)
}