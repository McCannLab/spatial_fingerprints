#' Get data formatted
#'
#' @param scale a logical. Should data se scaled?
#'
#' @export

get_data_ready <- function(scale = TRUE) {
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
  # Scale
  if (scale)
    df_dat[, -1] <- apply(df_dat[, -1], 2, scale)
  #
  df_dat <- data.frame(
        region = keepLetters(df_dat$id, 1:2), df_dat,
        stringsAsFactors = FALSE)
  #
  rownames(df_dat) <- NULL
  df_dat
}
