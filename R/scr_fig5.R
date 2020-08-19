#' Figures 5, 6 and 7
#'
#' Code to reproduce figures 5, 6 and 7 as well as Figure S8-S10.
#'
#' @param file path to results for individual performances as well as performances of all pairs and triplets.
#' @param meth a character string that describes the methods, used in the name of files that include the tables exported.
#' @param ind indices of the figures generated.

#' @export


scr_fig5 <- function(file = "output/res_lda_nb/all123/res_lda_123.rds", meth = "lda",
  ind = 5:7) {

  tmp <- lapply(readRDS(file), function(x) apply(x$mean, 3, function(y) mean(diag(y))))

  hh <- get_data_ready()
  res_ind <- data.frame(perf_ind = tmp[[1]], var = 0, var_wtn = 0, var_btw = 0,
    mean_distb = 0  # distance between centroids
)


  # intertia inter / intra + distance between points
  for (i in 3:19) {
    res_ind$var[i - 2] <- sum(hh[, i]^2)
    val1 <- hh[1:30, i]
    val2 <- hh[31:60, i]
    val3 <- hh[61:90, i]
    res_ind$var_wtn[i - 2] <- myvar(val1) + myvar(val2) + myvar(val3)
    m1 <- mean(hh[1:30, i])
    m2 <- mean(hh[31:60, i])
    m3 <- mean(hh[61:90, i])
    res_ind$var_btw[i - 2] <- 30 * (m1^2 + m2^2 + m3^2)
    res_ind$mean_distb[i - 2] <- mean((m1 - m2)^2, (m1 - m3)^2, (m2 - m3)^2)
    res_ind$dist2[i - 2] <- 1/90 * sum((val1 - m2)^2 + (val1 - m3)^2, (val2 -
      m1)^2 + (val2 - m3)^2, (val3 - m1)^2 + (val3 - m2)^2)
  }

  # Same work for pair of bio-tracers
  comb_bio2 <- data.frame(t(combn(3:19, 2)), perf = tmp[[2]])
  comb_bio2$max_perf1 <- apply(comb_bio2, 1, function(x) getmax(x[1:2], perf = res_ind$perf_ind))
  comb_bio2$mean_perf2 <- apply(comb_bio2, 1, function(x) getsum(x[1:2], perf = res_ind$perf_ind))
  comb_bio2$categ <- apply(comb_bio2, 1, function(x) get_categ2(x[1], x[2]))
  comb_bio2$mean_distb <- apply(comb_bio2, 1, function(x) getdist(hh[x[1:2]]))
  comb_bio2$dist2 <- apply(comb_bio2, 1, function(x) getdist2(hh[x[1:2]]))

  comb_bio2$max_dist <- apply(comb_bio2, 1, function(x) getmax(x[1:2], perf = res_ind$dist2))
  comb_bio2 <- cbind(comb_bio2, t(apply(comb_bio2, 1, function(x) inertia(hh[x[1:2]]))))

  names(comb_bio2)[1:2] <- paste0("bio", 1:2)
  comb_bio2$intersect <- apply(comb_bio2, 1, function(x) getintersect(hh[x[1:2]]))


  getval <- function(ids3, FUN = max) {
    FUN(comb_bio2$perf[apply(comb_bio2[, 1:2], 1, function(x) sum(x %in% ids3) ==
      2)])
  }

  comb_bio3 <- data.frame(t(combn(3:19, 3)), perf = tmp[[3]])
  names(comb_bio3)[1:3] <- paste0("bio", 1:3)
  comb_bio3$max_perf2 <- apply(comb_bio3, 1, function(x) getval(x[1:3]))
  comb_bio3$mean_perf2 <- apply(comb_bio3, 1, function(x) getval(x[1:3], FUN = mean))
  comb_bio3$mean_distb <- apply(comb_bio3, 1, function(x) getdist(hh[x[1:3]]))
  comb_bio3 <- cbind(comb_bio3, t(apply(comb_bio3, 1, function(x) inertia(hh[x[1:3]]))))
  comb_bio3$dist2 <- apply(comb_bio3, 1, function(x) getdist2(hh[x[1:2]]))
  comb_bio3$max_dist2 <- apply(comb_bio3, 1, function(x) getmax(x[1:3], perf = comb_bio2$dist2))
  comb_bio3$intersect <- apply(comb_bio3, 1, function(x) getintersect(hh[x[1:3]]))
  comb_bio3$max_intersect <- apply(comb_bio3, 1, function(x) getmax(x[1:3], perf = comb_bio2$intersect))
  # isoscapes names
  nmb <- gsub("n", "n-", gsub("_", ":", gsub("FA_", "", names(hh))))
  #
  tbl2 <- comb_bio2[rev(order(comb_bio2$perf))[1:10], 1:3]
  tbl2[, 1] <- nmb[tbl2[, 1]]
  tbl2[, 2] <- nmb[tbl2[, 2]]
  tbl2[, 3] <- round(tbl2[, 3], 3)
  tbl2 <- cbind(rank = 1:10, tbl2)
  #
  tbl3 <- comb_bio3[rev(order(comb_bio3$perf))[1:10], 1:4]
  for (i in 1:3) tbl3[, i] <- nmb[tbl3[, i]]
  tbl3[, 4] <- round(tbl3[, 4], 3)
  tbl3 <- cbind(rank = 1:10, tbl3)
  #
  output_dir()
  cat(kable(tbl2, row.names = FALSE), file = paste0("output/pairs10_", meth, ".md"), sep = "\n")
  cat(kable(tbl3, row.names = FALSE), file = paste0("output/triplets10_", meth, ".md"), sep = "\n")


  colr <- "#f63267"

  add_ticks <- function(at) axis(1, at = at, lwd = 0, lwd.ticks = 0.5, tck = -0.025,
    labels = NA)


  png("output/figs/fig5.png", width = 130, height = 60, units = "mm", res = 600)
  layout(matrix(1:3, ncol = 3), widths = c(1, 1, 0.35))
  par(las = 1, mar = c(3.2, 3.8, 1, 0), mgp = c(2.5, 0.7, 0), cex.axis = 0.8, cex.lab = 0.9)

  vc_pos <- rep(3, 17)
  vc_pos[c(8, 13, 9)] <- 1
  vc_pos[c(3, 14, 16)] <- 2
  vc_pos[c(6, 11, 15, 13)] <- 4
  rgy <- c(0.3, 0.6)
  # plot(apply(abs(res_ind[,-1]), 1, mean), res_ind[,1], pch = 19)
  palg <- c(rep("grey75", 3), rep("grey25", 14))
  dpalg <- sapply(palg, darken, 25)
  #
  plot(100 * res_ind$var_btw/89, res_ind$perf_ind, bg = palg, col = dpalg, pch = 21,
    cex = 0.9, ylim = rgy, xlim = 100 * c(0, 0.55), xlab = "", ylab = "Overall performance")
  par(mgp = c(2, 0.7, 0))
  title(xlab = "Inter-regions variance (%)")
  text(100 * res_ind$var_btw/89, res_ind$perf_ind, 1:17, pos = vc_pos, offset = 0.3,
    col = palg, cex = 0.8)
  add_ticks(seq(5, 55, 5))
  mtext("a", 3, at = 0, font = 2, cex = 0.8)
  ##
  vc_pos[c(3, 6, 10)] <- 1
  vc_pos[c(1)] <- 4
  vc_pos[c(8:9, 16)] <- 3
  plot(res_ind$mean_distb, res_ind$perf_ind, bg = palg, col = dpalg, pch = 21,
    ylim = rgy, xlab = "", ylab = "", cex = 0.9)
  title(xlab = "Mean distance between centroids")
  text(res_ind$mean_distb, res_ind$perf_ind, 1:17, pos = vc_pos, offset = 0.3,
    col = palg, cex = 0.8)
  add_ticks(seq(0.25, 3.25, 0.5))
  mtext("b", 3, at = 0, font = 2, cex = 0.8)

  par(mar = c(4.6, 0.4, 1, 0.8))
  plot(c(0, 1), c(1, 17), type = "n", axes = FALSE, ann = FALSE)
  for (i in seq_len(17)) {
    text(0.44, 18 - i, paste0(i, ". "), col = dpalg[i], pos = 2, cex = 0.8)
    text(0.2, 18 - i, nmb[i + 2], col = dpalg[i], pos = 4, cex = 0.8)
  }

  dev.off()

  msgSuccess_fig("5", "output/figs")





  add_ylab <- function(ylab) {
    par(mgp = c(2.8, 0.7, 0))
    title(ylab = ylab)
    par(mgp = c(1.8, 0.7, 0))
  }

  png("output/figs/fig6.png", width = 89, height = 100, units = "mm", res = 600)

  par(mfrow = c(2, 2), yaxs = "i", las = 1, mar = c(3.8, 4, 1.6, 0.5), mgp = c(1.8,
    0.7, 0), cex = 0.45)

  ## P1
  plot(comb_bio2$max_perf1, comb_bio2$perf, pch = 20, col = "grey10", xlab = "Best individual performance",
    ylab = "", xlim = c(0.3, 0.72), ylim = c(0.3, 0.85), cex = 0.6)
  add_ylab("Overall performance of a pair of bio-tracers")
  abline(a = 0, b = 1, lty = 3, col = colr, lwd = 1)
  mtext("a", 3, at = 0.3, font = 2, cex = 0.7)

  ## P2
  plot(comb_bio3$max_perf2, comb_bio3$perf, pch = 20, col = "grey10", xlab = "Best pair performance",
    ylab = "", xlim = c(0.3, 0.72), ylim = c(0.3, 0.85), cex = 0.6)
  add_ylab("Overall performance of a triplet")
  abline(a = 0, b = 1, lty = 3, col = colr, lwd = 1)
  mtext("b", 3, at = 0.3, font = 2, cex = 0.7)

  ## P3
  par(mar = c(4, 4, 1.4, 0.5))
  plot(0.5 * comb_bio2$mean_perf2, comb_bio2$perf, pch = 20, col = "grey10", xlab = c("Average overall performance",
    "of the bio-tracers in the pair"), ylab = "", xlim = c(0.3, 0.72), ylim = c(0.3,
    0.85), cex = 0.6)
  add_ylab("Performance of a pair of bio-tracers")
  abline(a = 0, b = 1, lty = 3, col = colr, lwd = 1)
  mtext("c", 3, at = 0.3, font = 2, cex = 0.7)

  ## P4
  plot(comb_bio3$mean_perf2, comb_bio3$perf, pch = 20, col = "grey10", xlab = c("Average over performance",
    "of the pair of bio-tracers in the triplet"), ylab = "", xlim = c(0.3, 0.72),
    ylim = c(0.3, 0.85), cex = 0.6)
  add_ylab("Performance of a triplet")
  abline(a = 0, b = 1, lty = 3, col = colr, lwd = 1)
  mtext("d", 3, at = 0.3, font = 2, cex = 0.7)

  dev.off()

  msgSuccess_fig("6", "output/figs")




  add_rsq <- function(rsq, x = 7.2, y = 0.32) {
    text(x, y, cex = 0.9, pos = 2, labels = expression(R^2 == ""), offset = 0)
    text(x, y, cex = 0.9, pos = 4, labels = paste0(format(100 * rsq, digit = 3),
      "%"), offset = 0.1)
  }

  png("output/figs/fig7.png", width = 89, height = 100, units = "mm", res = 600)

  layout(rbind(c(1, 2), 3))
  par(yaxs = "i", las = 1, mar = c(4, 4, 1.5, 0.5), mgp = c(2.5, 0.65, 0), cex = 0.45)

  prop2 <- comb_bio2$inert_btw/comb_bio2$inert_tot
  plot(100 * prop2, comb_bio2$perf, pch = 20, col = "grey10", xlab = "", ylab = "Overall performance of a pair of bio-tracer",
    ylim = c(0.3, 0.82), cex = 0.8)
  f <- fitexp2(comb_bio2, 100 * prop2, lty = 2, col = colr, lwd = 1)
  rsq <- 1 - deviance(f)/deviance(lm(perf ~ 1, data = comb_bio2))
  add_rsq(rsq, x = 45)
  title(xlab = "Inter-regions variance (%)")
  mtext("a", 3, at = 2, font = 2, cex = 0.7)
  add_ticks(seq(5, 50, 5))

  par(mgp = c(2.6, 0.65, 0))
  prop3 <- comb_bio3$inert_btw/comb_bio3$inert_tot
  plot(100 * prop3, comb_bio3$perf, pch = 20, col = "grey10", xlab = "", ylab = "Overall performance of a triplet",
    ylim = c(0.3, 0.82), cex = 0.8)
  f <- fitexp2(comb_bio3, 100 * prop3, lty = 2, col = colr, lwd = 1)
  rsq <- 1 - deviance(f)/deviance(lm(perf ~ 1, data = comb_bio3))
  add_rsq(rsq, x = 45)
  title(xlab = "Inter-regions variance (%)")
  mtext("b", 3, at = 2, font = 2, cex = 0.7)
  add_ticks(seq(5, 50, 5))


  plot(log10(comb_bio2$intersect), comb_bio2$perf, xlim = c(-2.2, -0.2), ylim = c(0.3,
    0.82), pch = 1, cex = 0.9, ylab = "Overall performance", xlab = "log10(region overlap)",
    lwd = 0.6)
  points(log10(comb_bio3$intersect), comb_bio3$perf, col = 1, pch = 20, cex = 0.7)
  legend("topright", legend = c("pair", "triplet"), pch = c(1, 19), bty = "n",
    pt.cex = c(1.2, 1), pt.lwd = 0.7)
  mtext("c", 3, at = -2.25, font = 2, cex = 0.7)
  f <- fitexp3(c(log10(comb_bio2$intersect), log10(comb_bio3$intersect)), c(comb_bio2$perf,
    comb_bio3$perf), lty = 2, col = colr, lwd = 1)
  rsq <- 1 - deviance(f)/deviance(lm(perf ~ 1, data = comb_bio3))
  add_rsq(rsq, x = -2.1)
  add_ticks(setdiff(seq(-2.2, -0.2, 0.1), c(-2, 0, 0.5)))

  dev.off()

  msgSuccess_fig("7", "output/figs")
  invisible(0)

}
