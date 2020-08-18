#' Figure S6
#'
#' Code to reproduce figure S6.
#'
#' @export


scr_figS6 <- function() {

  output_dir()
  msgInfo("Creating figure S6")

  png("output/figs/figS6.png", width = 10, height = 4.5, units = "in", res = 300)

  tmp <- read.table("output/res_ml/res_train_si.txt")
  names(tmp) <- c("nbb", "noise", "augment", "id_sample", "id_reg_test", "id_reg_true",
    "prob")
  res_ml <- aggregate(prob ~ noise * augment, mean, data = tmp[tmp$id_reg_test ==
    tmp$id_reg_true, ])
  un <- unique(res_ml$noise)
  pal <- colorRampPalette(graphicsutils::gpuPalette("cisl"))(length(un))
  pal[which(un == 0.01)] <- "#f63267"
  vc_lwd <- rep(1.6, length(un))
  vc_lwd[which(un == 0.01)] <- 3

  # Stable isotopes only
  par(mfrow = c(1, 3), las = 2, mar = c(4.2, 4, 2, 0.5), lend = 1)
  plot_figS6(res_ml, un, pal, vc_lwd, ylab = "Overall performance")
  legend("topleft", legend = unique(res_ml$noise), col = pal, lwd = vc_lwd, bty = "n", seg.len = 3)

  # Fatty acids only
  tmp <- read.table("output/res_ml/res_train_fa.txt")
  names(tmp) <- c("nbb", "noise", "augment", "id_sample", "id_reg_test", "id_reg_true",
    "prob")
  res_ml <- aggregate(prob ~ noise * augment, mean, data = tmp[tmp$id_reg_test ==
    tmp$id_reg_true, ])
  plot_figS6(res_ml, un, pal, vc_lwd, let = "(b)")

  # All bio-tracers
  tmp <- read.table("output/res_ml/res_train_all.txt")
  names(tmp) <- c("nbb", "noise", "augment", "id_sample", "id_reg_test", "id_reg_true",
    "prob")
  res_ml <- aggregate(prob ~ noise * augment, mean, data = tmp[tmp$id_reg_test ==
    tmp$id_reg_true, ])
  plot_figS6(res_ml, un, pal, vc_lwd, let = "(c)")

  dev.off()
  msgSuccess_fig("S6", "output/figs")

  invisible(0)
}



plot_figS6 <- function(res_ml, un, pal, vc_lwd, let = "(a)", ylab = "") {
  plot(c(0, 3.7), c(0.25, 1), type = "n", xlab = "Augmentation", ylab = ylab, axes = FALSE)
  axis(2, lwd = 0, lwd.ticks = 1)
  axis(1, at = log10(res_ml$augment), labels = res_ml$augment, lwd = 0, lwd.ticks = 1)
  for (i in seq_along(un)) {
    tmp2 <- res_ml[res_ml$noise == un[i], ]
    lines(log10(tmp2$augment), tmp2$prob, col = pal[i], lwd = vc_lwd[i])
  }
  abline(v = 3, col = "#f63267", lty = 2, lwd = 2)
  mtext(let, 3, at = 0, line = 0.3, cex = 1.2, font = 2, las = 1)
  box()
}
