#' Format data
#'
#' This is the function that has been used to format the raw data `df_all`.
#'
#' @param center a logical. Should data be centered?
#' @param scale a logical. Should data be scaled?
#' @param pca a logical. Should column be transformed and order via a PCA. Note that in this case, data are scaled.
#'
#' @export

get_data_ready <- function(center = TRUE, scale = TRUE, pca = FALSE) {
  df_dat <- spatialfingerprints::df_all
  ## storage not needed as id is explicit
  # df_dat$Storage <- c("Lif", "Frozen")[2 - grepl("L$", df_dat$id)]
  ## rename
  names(df_dat) <- gsub("cal_", "", names(df_dat))
  ## subset column
  df_dat <- df_dat[-c(3, 5:6, 8)]
  ## keep only frozen
  # df_dat <- df_dat[df_dat$Storage == "Frozen", ]
  df_dat <- df_dat[!grepl("L$", df_dat$id), ]
  #
  if (pca) {
    df_dat[, -1] <- prcomp(df_dat[, -1], center = center, scale. = scale)$x
  } else {
    # Scale
    df_dat[, -1] <- apply(df_dat[, -1], 2, scale, center, scale)
  }
  #
  df_dat <- data.frame(
        region = keepLetters(df_dat$id, 1:2), df_dat,
        stringsAsFactors = FALSE)
  #
  rownames(df_dat) <- NULL
  df_dat
}
