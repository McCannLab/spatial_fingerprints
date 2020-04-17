#' Core functions to asses overall performance of LDA and KS.
#'
#' @param ndistr number of sample used to build baseline biotracer distributions.
#' @param nsample number of
#' @param noise noise to be added.
#' @param col_ids ids (name or integer) of columns to be used (select biotracers).
#' @param df_dat biotracer data sets. Note that by default, it loads `df_all()`
#' (see [get_data_ready()]).
#' @param toprob should the output be a probability.
#'
#' @details
#' `ndistr + nsample` cannot exceed 30 with data in `df_all`
#'
#' @export

## Use lda (MASS) https://en.wikipedia.org/wiki/Linear_discriminant_analysis
scr_lda <- function(ndistr = 20, nsample = 10, noise = 0,
  col_ids = 3:19, df_dat = NULL, toprob = TRUE) {
  if (is.null(df_dat))
    df_dat <- get_data_ready()
  #
  ls_val <- split(seq_along(df_dat$region), df_dat$region)
  ngeo <- length(unique(df_dat$region))
  #
  id <- lapply(ls_val, sample, size = ndistr + nsample)
  sq <- seq_len(nsample)
  ids <- lapply(id, `[`, sq)
  idd <- unlist(lapply(id, `[`, -sq))
  #
  nms <- c("region", names(df_dat)[col_ids])
  tmp <- df_dat[nms]
  train <- tmp[idd, ]
  # add noise to distribution for the training
  if (noise) {
    for (i in nms[-1])
      train[i] <- train[i] + rnorm(ndistr, sd = noise)
  }
  fml <- paste0("region ~ ", paste(nms[-1], collapse = " + "))
  mod <- lda(as.formula(fml), prior = rep(1, ngeo)/3, data = train)
  ##
  prd <- predict(mod, df_dat[unlist(ids), nms[-1], drop = FALSE])$posterior
  ## code below combine predictions of several samples (i.e. a couple of fish)
  ## matching sample/country work when the rownames have been re-numbered from
  ## 1 to nrow(df_dat)
  out <- do.call(rbind, lapply(lapply(lapply(ids, as.character), function(x) prd[x, ,drop = FALSE]), apply, 2, function(x) -sum(log(x))))
  ##
  if (toprob)
    toprob(out) else out
}

#' @describeIn scr_lda core function for KS
#' @export
scr_ks <- function(ndistr = 25, nsample = 5, noise = 0,
  col_ids = c("d13C", "d15N", "d34S"), df_dat = NULL, toprob = TRUE) {
  if (is.null(df_dat))
    df_dat <- get_data_ready()

  ls_val <- split(seq_along(df_dat$country), df_dat$country)
  ngeo <- length(unique(df_dat$country))
  ##
  tmp <- df_dat[col_ids]
  id <- lapply(ls_val, sample, size = ndistr + nsample)
  idd <- seq_len(ndistr)
  # use the same bandwith for all data
  H <- Hpi(tmp)
  out <- matrix(0, ngeo, ngeo)
  ## over the 3 regions
  for (j in seq_len(ngeo)) {
    train <- tmp[id[[j]][idd], ]
    if (noise)
      train <- add_noise(train, sd = noise)
    for (l in seq_len(ngeo)) {
      smpl <- tmp[id[[l]][-idd], ]
      ## use LL NB: for 1 sample (1 row) test its origin, i.e. proba it comes from
      ## regions (all columns)
      out[l, j] <- -sum(log(kde(train, H, eval.points = smpl)$estimate))
    }
  }
  if (toprob)
    toprob(out) else out
}



# convert a matrix of -log likelihood into a proba for each state
# NB: FOR EACH ELEMENR row i column j
# $exp(-m_{i,j})/\sum_j{exp(-m_{i,j})}$
##
toprob <- function(x) {
  tmp <- exp(-x)
  # NB: call to t() required, otherwise results are in the wrong order
  t(apply(tmp, 1, function(x) x/sum(x)))
}



# add noise
add_noise <- function(train, mean = 0, sd = 1) {
  apply(train, 2, function(x) x + rnorm(length(x), mean, sd))
}

add_noise_v <- function(train, mean = 0, sd = 1) {
  train + rnorm(length(train), mean, sd)
}

scr_ks_ind <- function(ndistr = 20, nsample = 5, noise = 0, col_ids = 3:19, df_dat = NULL, toprob = TRUE) {
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

