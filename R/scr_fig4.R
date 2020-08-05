#' Figure 4
#'
#' Code to reproduce figure 4.
#'
#' @param path Path to results files for LDA and NBC.
#' @param file_ml Path to results file (a rds file) for MLP.
#' @param files_pca A vector of 3 rds files (1 per method) corresponding to the results obtained once a PCA was applied to the data set.
#'
#' @export

scr_fig4 <- function(path = "output/res_lda_nb/nbio/non_pca", file_ml = "output/res_f/res_combn_ml_nbio.rds",
  files_pca = c("output/res_f/nbio_pca_lda_01.rds", "output/res_f/nbio_pca_nb_01.rds", "output/res_f/ml_nbio_pca.rds")) {

  ## Helper functions
  get_res_bb <- function(file) {
    raw <- readRDS(file)
    lapply(raw, function(x) apply(x$mean, 3, function(x) mean(diag(x))))
  }

  get_res_pca <- function(file) {
    raw <- readRDS(file)
    lapply(raw, function(x) apply(x, c(1, 2), mean))
  }

  get_res_sd <- function(file) {
    raw <- readRDS(file)
    lapply(raw, function(x) apply(x$sd, c(1, 2), mean))
  }

  addlet <- function(let, x = 1) mtext(let, 3, at = x, cex = 1, font = 2)

  addaxes <- function() {
    axis(1, at = seq(2, 16, 2), labels = NA, lwd = 0, lwd.ticks = 0.25, tck = -0.025)
    axis(1, at = seq(1, 17, 2), labels = seq(1, 17, 2), lwd = 0, lwd.ticks = 0.5)
    axis(2, lwd = 0, lwd.ticks = 0.5)
    box()
  }

  ## Read files
  files_lda <- sprintf(paste0(path, "/nbio_lda_%02d.rds"), 1)
  files_nb <- sprintf(paste0(path, "/nbio_lda_%02d.rds"), 1)
  pca_lda <- unlist(lapply(get_res_pca(files_pca[1L]),
    function(x) mean(diag(x))))
  pca_nb <- unlist(lapply(get_res_pca(files_pca[2L]),
    function(x) mean(diag(x))))
  # MLP
  tmp <- readRDS(files_pca[3L])
  res_ml <- tmp[tmp$id_reg_test == tmp$id_reg_true, ]
  ml_pca <- aggregate(prob ~ nbio, mean, data = res_ml)
  tmp2 <- readRDS("output/res_f/res_ml_nbio.rds")
  res_ml2 <- tmp2[tmp2$id_reg_test == tmp2$id_reg_true, ]
  ml_reg <- aggregate(prob ~ nbio * id_comb, mean, data = res_ml2)
  pca_ml <- split(ml_reg$prob, ml_reg$nbio)

  sq_bt <- seq_len(17)

  output_dir("output/figs")
  msgInfo("Creating figure 4")
  png("output/figs/fig4.png", width = 183, height = 72, units = "mm", res = 600)

  par(mfrow = c(1, 3), las = 1, mar = c(4, 3.2, 1.5, 0.4),
    mgp = c(2.25, 0.6, 0))
  plot(range(sq_bt), c(0.33, 1), type = "n", axes = FALSE, xlab = "",
    ylab = "Overall performance")
  boxplot(get_res_bb(files_lda[[1L]]), col = "grey95", add = TRUE, pch = 19, border = "grey55",
    lwd = 0.8, cex = 0.5, axes = FALSE)
  points(sq_bt, pca_lda, col = 1, pch = 19, cex = 1)
  lines(sq_bt, pca_lda, col = 1, lwd = 0.7)
  addaxes()
  addlet("a")

  plot(range(sq_bt), c(0.33, 1), type = "n", xlab = "Number of bio-tracers combined",
    ylab = "", axes = FALSE)
  boxplot(get_res_bb(files_nb[[1L]]), col = "grey95", add = TRUE, xes = FALSE,
    pch = 19, border = "grey55", lwd = 0.8, cex = 0.5)
  points(sq_bt, pca_nb, col = 1, pch = 19, cex = 1)
  lines(sq_bt, pca_nb, col = 1, lwd = 0.7)
  addaxes()
  addlet("b")

  plot(range(sq_bt), c(0.33, 1), type = "n", xlab = "", ylab = "", axes = FALSE)
  boxplot(pca_ml, col = "grey95", add = TRUE, pch = 19, border = "grey55", lwd = 0.8,
    cex = 0.5, axes = FALSE)
  points(sq_bt, ml_pca[, 2], col = 1, pch = 19, cex = 1)
  lines(sq_bt, ml_pca[, 2], col = 1, lwd = 0.7)
  addaxes()
  addlet("c")

  dev.off()

  msgSuccess_fig("4", "output/figs")
  invisible(0)

}

