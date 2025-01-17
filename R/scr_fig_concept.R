#' Figure 1
#'
#' Code to reproduce figure 1
#'
#' @details
#' Note that the Salmon silhouette was added with Inkscape and is available at
#' <http://phylopic.org/image/3c098bb8-4158-4777-b567-80e48049409c/>.
#'
#' @export

scr_fig_concept <- function() {

  msgInfo("Creating figure 1")

  mp <- orthomap(centre = c(45, 175))
  dev.off()
  pal <- gpuPalette("insileco")[2:4]
  mp2 <- st_as_sf(mp)

  xpt <- c(0.445, 0.131, -0.284)
  ypt <- c(0.244, 0.225, 0.157)
  xc <- .1
  yc <- -.15
  muC <- c(3.5, 4, 6.5)
  siC <- c(0.7, 0.9, 1.8)
  muN <- c(2, 5, 6)
  siN <- c(0.8, 1, 1.6)
  sq <- seq(0, 10, 0.05)
  mat <- matrix(c(1, 1, 2, 4, 0, 3), 2)


  output_dir()
  png("output/fig1.png", width = 8, height = 5, units = "in", res = 600)

  layout(mat, widths = c(.9, 1, .3), heights = c(.3, 1))
  ## MAP
  par(mar = c(0, 0, 0, 0), cex.axis = 1.2, cex.lab = 1.6,
    bg = "transparent", fg = "gray20", col = "gray20",  col.lab = "gray20",
    col.axis = "gray20", lwd = 2)
  plot(st_geometry(mp2), col = "grey80", border = "grey40", lwd = .3)
  for (i in seq_along(xpt)) lines(c(xpt[i], xc), c(ypt[i], yc))
  points(c(0.445, 0.131, -0.284), c(0.244, 0.225, 0.157), cex = 3, pch = 19, col = pal)
  text(xc, yc, labels = "?", pos = 1, cex = 4, offset = .8)

  par(bty = "l", las = 1, lend = 1, xaxs = "i", yaxs = "i")
  ##
  xbio <- c(0, 10)
  ybio <- c(-.03, 0.6)
  par(mar = c(0, 4, .6, 0))
  plot(xbio, ybio, type = "n", axes = FALSE, ann = FALSE)
  abline(v = 3.8, lwd = .9, lty = 2)
  for (i in 1:3) lines(sq, dnorm(sq, muC[i], siC[i]), col = pal[i], lwd = 2)
  #
  par(mar = c(4, 0, 0, .6))
  plot(ybio, xbio, type = "n", axes = FALSE, ann = FALSE)
  abline(h = 5.6, lwd = .9, lty = 2)
  for (i in 1:3) lines(dnorm(sq, muN[i], siN[i]), sq, col = pal[i], lwd = 2)

  ##
  par(mar = c(4, 4, 0, 0), bty = "l", mgp = c(2.1, .6, 0))
  sigm <- matrix(0, 2, 2)
  for (i in 1:3) {
    diag(sigm) <- c(siC[i], siN[i])
    val <- matrix(dmvnorm(expand.grid(sq, sq), c(muC[i], muN[i]),
      sigm, log = FALSE), length(sq))
    val2 <- scaleWithin(val, 256, 10^-6, 10^-1)
    pal2 <- paste0(pal[i], sprintf("%02X", 0:255))
    if (i == 1) {
      image(sq, sq, val2, col = pal2, xlab = "Biotracer 1", ylab = "Biotracer 2")
    } else image(sq, sq, val2, col = pal2, add = TRUE)
  }
  abline(v = 3.8, h = 5.6, lwd = .9, lty = 2)

  dev.off()

}






other_version <- function() {
  mp <- orthomap(centre = c(45, 175))
  dev.off()
  mp2 <- st_as_sf(mp)

  png("map2d.png", width = 100, height = 100, units = "mm", res = 600)
  # svg(filename, width = 7, height = 4, pointsize = 12)
  pal <- c("#ee2485", 1, "#3fb3b2")

  ## MAP
  par(mar = c(0, 0, 0, 0), cex.axis = 1.2, cex.lab = 1.6,  bg = "transparent", fg =
  "gray20", col = "gray20",  col.lab = "gray20", col.axis = "gray20",
  lwd = 2)
  plot(st_geometry(mp2), col = "#ebc5d7", border = "#ee2485", lwd = .3)
  xpt <- c(-0.1236271, 0.1844303, 0.5037306)
  ypt <- c(0.1290981, 0.2393756, 0.2843869)
  points(xpt, ypt, cex = 2, pch = 21, bg = c("#ebc5d7", "grey50", "#79b1b1"), col = pal, lwd = 3)

  # for (i in seq_along(xpt)) lines(c(xpt[i], xc), c(ypt[i], yc))
  # points(c(0.445, 0.131, -0.284), c(0.244, 0.225, 0.157), cex = 3, pch = 19, col = pal)
  # text(xc, yc, labels = "?", pos = 1, cex = 4, offset = .8)
  ##
  dev.off()
}

# Matrices for illustration purposes
fig_part3 <- function() {

  mat1 <- apply(replicate(2000, find_origin("nb", ndistr = 20, nsample = 2,
    col_ids =  13)), c(1, 2), mean)
  mat2 <- apply(replicate(2000, find_origin("nb", ndistr = 20, nsample = 2,
    col_ids =  7)), c(1, 2), mean)
  mat3 <- apply(replicate(2000, find_origin("nb", ndistr = 22, nsample = 6,
    col_ids =  c(13, 7))), c(1, 2), mean)
  pal <- colorRampPalette(c("white", "black"))(100)
  png("output/fig1c.png", width = 3, height = 5, units = "in", res = 600)
  par(mfrow = c(3, 1), mar = c(2, 1, 1, 1), bg = "transparent")
  graphicsutils::image2(scaleWithin(mat1, 100, 0, 1), color_scale = pal)
  mat_text(format(mat1, digits = 2))
  graphicsutils::image2(scaleWithin(mat2, 100, 0, 1), color_scale = pal)
  mat_text(format(mat2, digits = 2))
  graphicsutils::image2(scaleWithin(mat3, 100, 0, 1), color_scale = pal)
  mat_text(format(mat3, digits = 2))
  dev.off()
  c(mean(diag(mat1)), mean(diag(mat2)), mean(diag(mat3)))
}


mat_text <- function(mat, ...) {
  seqx <- seq(0, 1, length.out = NROW(mat))
  seqy <- seq(0, 1, length.out = NCOL(mat))
  for (i in seq_len(NROW(mat))) {
    for (j in seq_len(NCOL(mat))) {
      col <- ifelse(as.numeric(mat[i, j]) < .5, "black", "white")
      text(seqx[i], 1 - seqy[j], labels = mat[i, j], col = col, cex = 1.1, ...)
    }
  }
}