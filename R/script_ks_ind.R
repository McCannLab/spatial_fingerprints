script_ks_ind <- function(ndistr = 20, nsample = 5, noise = 0, col_ids = 3:19, df_dat = NULL, toprob = TRUE) {
  if (is.null(df_dat))
    df_dat <- get_data_ready()

  ls_val <- split(seq_along(df_dat$region), df_dat$region)
  ngeo <- length(unique(df_dat$region))
  ##
  tmp <- df_dat[col_ids]
  id <- lapply(ls_val, sample, size = ndistr + nsample)
  idd <- seq_len(ndistr)
  # use the same bandwith for all data
  out <- matrix(0, ngeo, ngeo)
  for (k in seq_along(col_ids)) {
    for (j in seq_len(ngeo)) {
      tmp2 <- tmp[id[[j]][idd], k]
      if (noise)
        tmp2 <- add_noise_v(tmp2, sd = noise)
      train <- density(tmp2)
      for (l in seq_len(ngeo)) {
        smpl <- tmp[id[[l]][-idd], k]
        ## use LL
        out[l, j] <- out[l, j] + isos_likelihood(smpl, train)
      }
    }
  }
  if (toprob)
    toprob(out) else out
}
