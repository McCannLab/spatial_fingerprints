#' Core functions to asses overall performance of LDA, NB and ML.
#'
#' @param method statistical approach, either "lda", "nb" or "ml".
#' @param df_dat biotracer data sets. Note that by default, it loads `df_all()`
#' (see [get_data_ready()]).
#' @param ndistr number of sample used to build baseline biotracer distributions.
#' @param nsample number of samples.
#' @param noise noise to be added.
#' @param col_ids ids (name or integer) of columns to be used (select bio-tracers).
#' @param toprob should the output be a probability.
#'
#' @details
#' * "lda" : linear discriminant analysis,
#' * "nb" : naive bayesian,
#' * "ml" : machine learning,
#' * `ndistr + nsample` cannot exceed 30 with data in `df_all`.
#'
#' @export
#' @examples
#' df_dat <- get_data_ready()
#' find_origin("nb", col_ids = 3:4, df_dat = df_dat)

find_origin <- function(method = c("lda", "nb", "ml"), df_dat, ndistr = 20,
  nsample = 10, noise = 0, col_ids = 3:19, toprob = TRUE) {

    method <- match.arg(method)
    ls_val <- split(seq_len(nrow(df_dat)), df_dat$region)
    ngeo <- length(ls_val)

    out <- switch(method,
      lda = scr_lda(df_dat, ndistr, nsample, noise, col_ids),
      nb = scr_nb(df_dat, ndistr, nsample, noise, col_ids),
      ml = stop("implemented in Julia")
    )

    if (toprob) toprob(out) else out
}


#' @describeIn find_origin core function for LDA.
#' @export
scr_lda <- function(df_dat, ndistr = 20, nsample = 10, noise = 0,
  col_ids = 3:19) {

  ## Use lda (MASS) <https://en.wikipedia.org/wiki/Linear_discriminant_analysis>
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
  out
}




#' @describeIn find_origin core function for the Naive Bayesian (NB) classifier.
#' @export
scr_nb <- function(df_dat, ndistr = 20, nsample = 5, noise = 0, col_ids = 3:19) {

  ls_val <- split(seq_along(df_dat$region), df_dat$region)
  ngeo <- length(unique(df_dat$region))
  ##
  tmp <- df_dat[col_ids]
  # randomize training set and test set
  id <- lapply(ls_val, sample, size = ndistr + nsample)
  idd <- seq_len(ndistr)
  # use the same bandwith for all data
  out <- matrix(0, ngeo, ngeo)
  for (k in seq_along(col_ids)) {
    for (j in seq_len(ngeo)) {
      # training set
      tmp2 <- tmp[id[[j]][idd], k]
      if (noise)
        tmp2 <- add_noise_v(tmp2, sd = noise)
      train <- density(tmp2, adjust = 2)
      for (l in seq_len(ngeo)) {
        smpl <- tmp[id[[l]][-idd], k]
        ## LL for all samples and for NB LL is added for all bio-tracers
        out[l, j] <- out[l, j] + isos_likelihood(smpl, train)
      }
    }
  }
  colnames(out) <- rownames(out) <- c("CA", "RS", "US")
  out
}

