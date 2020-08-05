#' Analysis pipeline
#'
#' Run all analyses.
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
  #

  # for (i in 1:20) {
  #   out <- simu_nbio("lda", mx_cb = 1000, nrep = 1000, nsample = i,
  #     ndistr = 20)
  #   saveRDS(out, file = paste0("lda_nbio_", nsample, ".rds")
  #   out <- simu_nbio("nb", mx_cb = 1000, nrep = 1000, nsample = i,
  #     ndistr = 20)
  #   saveRDS(out, file = paste0("lda_nbio_", nsample, ".rds")
  # }

  msgInfo("Main Figures...")
  scr_fig_concept()
  scr_fig2()
  scr_fig3()


  msgInfo("SI Figures...")
  scr_theory()
  #
  invisible(NULL)
}