#' Figure 3
#'
#' Code to reproduce figure 3 (noise addition).
#'
#' @param path Path to results files for LDA and NBC.
#' @param file_ml Path to results file (a rds file) for MLP.
#'
#' @export

scr_fig3 <- function(path = "output/res_lda_nb/noise/", file_ml = "output/res_f/res_combn_ml_nbio.rds") {

  ## Helper functions
  get_res <- function(file) {
    # read files
    raw <- readRDS(file)
    # remove empty elements
    raw <- Filter(Negate(is.null), raw)
    lapply(raw, function(x) apply(x$mean, c(1, 2), mean))
  }

  addlet <- function(let, x = 1) mtext(let, 3, at = x, cex = 1, font = 2)

  addaxesN <- function() {
    axis(1, at = setdiff(1:26, seq(1, 26, 5)), labels = NA, lwd = 0, lwd.ticks = 0.25,
      tck = -0.025)
    axis(1, at = seq(1, 26, 5), labels = c(TeX("$10^{-4}$"), TeX("$10^{-3}$"),
      TeX("$10^{-2}$"), TeX("$10^{-1}$"), TeX("$1$"), TeX("$10$")), lwd = 0,
      lwd.ticks = 0.5)
    axis(2, lwd = 0, lwd.ticks = 0.5)
    box()
  }

  ## reading files
  idf <- c(1:2, 5, 10, 15)
  nf <- length(idf)
  msgInfo("Reading files for figure 3")
  files_lda <- sprintf(paste0(path, "/noise_lda_%02d.rds"), idf)
  files_nb <- sprintf(paste0(path, "/noise_nb_%02d.rds"), idf)
  
  # the following steps required some RAM so there were commented out and 
  # an intermediate file has been saved. 
  # tmp <- readRDS("output/res_ml/ml_noise.rds")
  # get data for which region tested in the true one
  # res_ml <- tmp[tmp$id_reg_test == tmp$id_reg_true, ]
  # ml_reg <- aggregate(prob ~ nbio * noise, mean, data = res_ml)
  # saveRDS(ml_reg, "output/res_ml/ml_reg_noise.rds")
  ml_reg <- readRDS("output/res_ml/ml_reg_noise.rds")

  ## Figure
  output_dir("output/figs")
  msgInfo("Creating figure 3")
  sqn <- 1:26
  pal <- colorRampPalette(c("grey5", "grey75"))(nf)


  png("output/figs/fig3.png", width = 183, height = 70, units = "mm", res = 600)

  par(mfrow = c(1, 3), las = 1, mar = c(4, 3.2, 1, 0.4), mgp = c(2.25, 0.6, 0))

  ## LDA
  plot(range(sqn), c(0.33, 1), type = "n", axes = FALSE, xlab = "",
    ylab = "Overall performance")
  for (i in seq_len(nf)) {
    val <- unlist(lapply(get_res(files_lda[[i]]), function(x) mean(diag(x))))
    lines(sqn, val, pch = 19, col = pal[i], lwd = 0.7)
    points(sqn, val, pch = 19, col = pal[i], cex = 0.8)
  }
  addlet("a")
  addaxesN()
  # NBC
  plot(range(sqn), c(0.33, 1), type = "n", xlab = "Noise level", ylab = "", axes = FALSE)
  for (i in seq_len(nf)) {
    val <- unlist(lapply(get_res(files_nb[[i]]), function(x) mean(diag(x))))
    lines(sqn, val, pch = 19, col = pal[i], lwd = 0.7)
    points(sqn, val, pch = 19, col = pal[i], cex = 0.8)
  }
  addlet("b")
  addaxesN()
  ## MLP
  plot(range(sqn), c(0.33, 1), type = "n", xlab = "", ylab = "", axes = FALSE)
  # ml
  for (i in seq_along(idf)) {
    # cx <- ml_reg$noise[ml_reg$nbio == idf[i]]
    cx <- sqn
    cy <- ml_reg$prob[ml_reg$nbio == idf[i]]
    lines(cx, cy, pch = 19, col = pal[i], lwd = 0.7)
    points(cx, cy, pch = 19, col = pal[i], cex = 0.8)
  }
  addlet("c")
  addaxesN()
  #
  legend("topright", legend = idf, col = pal, pch = 19, bty = "n")

  dev.off()

  msgSuccess_fig("3", "output/figs")
  invisible(0)

}





#' @describeIn scr_fig3 Code to reproduce figure S7.
#' @export
scr_figS7 <- function() {

  ## Helper functions
  get_res <- function(file) {
    raw <- readRDS(file)
    # remove empty elements
    raw <- Filter(Negate(is.null), raw)
    lapply(raw, function(x) apply(x$mean, c(1, 2), mean))
  }

  addaxesD <- function() {
    sq <- seq(1, 21, 5)
    axis(1, at = setdiff(1:21, sq), labels = NA, lwd = 0, lwd.ticks = 0.25, tck = -0.025)
    axis(1, at = sq, labels = sq + 4, lwd = 0, lwd.ticks = 0.5)
    axis(2, lwd = 0, lwd.ticks = 0.5)
    box()
  }

  addlet <- function(let, x = 1) mtext(let, 3, at = x, cex = 1, font = 2)

  add_vl <- function(v = 16) abline(v = v, lwd = 1.2, col = "#f63267", lty = 2)

  ## Very similar => fig S
  idf <- c(1:2, 5, 10, 15)
  files_lda <- sprintf("output/res_lda_nb/ndistr/ndistr_lda_%02d.rds", idf)
  files_nb <- sprintf("output/res_lda_nb/ndistr/ndistr_nb_%02d.rds", idf)
  nf <- length(idf)
  tmp <- readRDS("output/res_ml/ml_ndistr.rds")
  res_ml <- tmp[tmp$id_reg_test == tmp$id_reg_true, ]
  ml_reg <- aggregate(prob~nbio*ndistr, mean, data = res_ml)


  # tmp <- readRDS('output/res_f/res_ml_nbio.rds') res_ml <- tmp[tmp$id_reg_test ==
  # tmp$id_reg_true, ]) ml_sam <- aggregate(prob~nbio, mean, data = res_ml)

  output_dir("output/figs")
  msgInfo("Creating figure S7")
  png("output/figs/figS7.png", width = 183, height = 70, units = "mm", res = 600)

  par(mfrow = c(1, 3), las = 1, mar = c(4, 3.2, 1, 0.4), mgp = c(2.25, 0.6, 0))

  sqd <- 1:21

  ## lda
  plot(range(sqd), c(0.33, 1), type = "n", xlab = "", ylab = "Overall performance",
    axes = FALSE)
  pal <- colorRampPalette(c("grey10", "grey70"))(nf)
  for (i in seq_len(nf)) {
    val <- unlist(lapply(get_res(files_lda[[i]]), function(x) mean(diag(x))))
    val <- val[-c(1, length(val))]
    lines(sqd, val, pch = 19, col = pal[i], lwd = 0.7)
    points(sqd, val, pch = 19, col = pal[i], cex = 0.8)
  }
  addlet("a")
  addaxesD()
  add_vl()

  # nb
  plot(range(sqd), c(0.33, 1), type = "n", xlab = "Size of the training set", ylab = "", axes = FALSE)
  for (i in seq_len(nf)) {
    val <- unlist(lapply(get_res(files_nb[[i]]), function(x) mean(diag(x))))
    val <- val[-c(1, length(val))]
    lines(sqd, val, pch = 19, col = pal[i], lwd = 0.7)
    points(sqd, val, pch = 19, col = pal[i], cex = 0.8)
  }
  addlet("b")
  addaxesD()
  add_vl()

  # ml
  plot(c(5, 25), c(0.33, 1), type = "n", xlab = "", ylab = "")
  for (i in seq_along(idf)) {
    cx <- ml_reg$ndistr[ml_reg$nbio == idf[i]]
    cy <- ml_reg$prob[ml_reg$nbio == idf[i]]
    lines(cx, cy, pch = 19, col = pal[i], lwd = 0.7)
    points(cx, cy, pch = 19, col = pal[i], cex = 0.8)
  }
  addlet("c")
  axis(1, at = setdiff(5:25, seq(5, 25, 5)), labels = NA, lwd = 0,
    lwd.ticks = .25)
  add_vl(v = 20)
  #
  legend("bottomright", legend = idf, col = pal, pch = 19, bty = "n",
    ncol = 5, cex = 1.025)

  dev.off()

  msgSuccess_fig("S7", "output/figs")
  invisible(0)

}