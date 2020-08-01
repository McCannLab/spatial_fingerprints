scr_fig4 <- function() {

get_res <- function(file) {
  raw <- readRDS(file)
  # remove empty elements
  raw <- Filter(Negate(is.null), raw)
  lapply(raw, function(x) apply(x$mean, c(1, 2), mean))
}

addlet <- function(let, x = 1) mtext(let, 3, at = x, cex = 1, font =2)


addaxesN <- function() {
  axis(1, at = setdiff(1:26, seq(1, 26, 5)), labels = NA, lwd = 0, lwd.ticks = .25, tck = -.025)
  axis(1, at = seq(1, 26, 5), labels = c(TeX("$10^{-4}$"), TeX("$10^{-3}$"), TeX("$10^{-2}$"), TeX("$10^{-1}$"), TeX("$1$"), TeX("$10$")), lwd = 0, lwd.ticks = .5)
  axis(2,  lwd = 0,  lwd.ticks = .5)
  box()
}



idf <- c(1:3, 5, 10, 15)
files_lda <- sprintf("output/res_lda_nb/noise/noise_lda_%02d.rds", idf)
files_nb <- sprintf("output/res_lda_nb/noise/noise_nb_%02d.rds", idf)
nf <- length(idf)
#
# tmp <- readRDS('output/res_f/res_ml_nbio.rds')
# res_ml <- tmp[tmp$id_reg_test == tmp$id_reg_true, ]
# ml_reg <- aggregate(prob~nbio*id_reg_true, mean, data = res_ml)
# ml_sam <- aggregate(prob~nbio, mean, data = res_ml)

png('output/figs/fig4.png',  width = 183, height = 70, units = "mm", res = 600)

par(mfrow = c(1, 3), las = 1, mar = c(4, 3.2, 1, .4), mgp = c(2.25, .6, 0))

sqn <- 1:26

## lda
plot(c(1, 26), c(0.33, 1), type = "n", xlab = "", ylab = "Overall performance",
axes = FALSE)
pal <- colorRampPalette(c('grey5', 'grey75'))(nf)
for (i in seq_len(nf)) {
  val <- unlist(lapply(get_res(files_lda[[i]]), function(x) mean(diag(x))))
  lines(sqn, val, pch = 19, col = pal[i], lwd = .7)
  points(sqn, val, pch = 19, col = pal[i], cex = .8)
}
addlet("a")
addaxesN()

# nb
plot(range(sqn), c(0.33, 1), type = "n", xlab = "Noise", ylab = "", axes = FALSE)
for (i in seq_len(nf)) {
  val <- unlist(lapply(get_res(files_nb[[i]]), function(x) mean(diag(x))))
  lines(sqn, val, pch = 19, col = pal[i], lwd = .7)
  points(sqn, val, pch = 19, col = pal[i], cex = .8)
}
addlet("b")
addaxesN()
# ml
plot(range(sqn), c(0.33, 1), type = "n", xlab = "", ylab = "")
addlet("c")
# points(ml_sam[, 1], ml_sam[, 2], pch = 19)
# lines(ml_sam[, 1], ml_sam[, 2], pch = 19)
# addlet("f")

#
legend("topright", legend = idf, col = pal, pch = 19, bty = "n")

dev.off()

}






scr_figSX <- function() {

  get_res <- function(file) {
    raw <- readRDS(file)
    # remove empty elements
    raw <- Filter(Negate(is.null), raw)
    lapply(raw, function(x) apply(x$mean, c(1, 2), mean))
  }

  addaxesD <- function() {
    sq <- seq(1, 23, 3)
    axis(1, at = setdiff(1:23, sq), labels = NA, lwd = 0, lwd.ticks = .25, tck = -.025)
    axis(1, at = sq, labels = sq + 3,lwd = 0, lwd.ticks = .5)
    axis(2,  lwd = 0,  lwd.ticks = .5)
    box()
  }

  addlet <- function(let) text(1, .98, let, cex = 1.2, font = 2)


## Very similar => fig SX


idf <- c(1:3, 5, 10, 15)
files_lda <- sprintf("output/res_lda_nb/ndistr/ndistr_lda_%02d.rds", idf)
files_nb <- sprintf("output/res_lda_nb/ndistr/ndistr_nb_%02d.rds", idf)
nf <- length(idf)


# tmp <- readRDS('output/res_f/res_ml_nbio.rds')
# res_ml <- tmp[tmp$id_reg_test == tmp$id_reg_true, ]
# ml_reg <- aggregate(prob~nbio*id_reg_true, mean, data = res_ml)
# ml_sam <- aggregate(prob~nbio, mean, data = res_ml)



png('output/figs/figSX.png', width = 183, height = 70, units = "mm", res = 600)

par(mfrow = c(1, 3), las = 1, mar = c(4, 3.2, 1, .4), mgp = c(2.25, .6, 0))

sqd <- 1:23

## lda
plot(range(sqd), c(0.33, 1), type = "n", xlab = "", ylab = "Overall performance",
axes = FALSE)
pal <- colorRampPalette(c('grey10', 'grey70'))(nf)
for (i in seq_len(nf)) {
  val <- unlist(lapply(get_res(files_lda[[i]]), function(x) mean(diag(x))))
  lines(sqd, val, pch = 19, col = pal[i], lwd = .7)
  points(sqd, val, pch = 19, col = pal[i], cex = .8)
}
addlet("a")
addaxesD()

# nb
plot(range(sqd), c(0.33, 1), type = "n", xlab = "Noise", ylab = "", axes = FALSE)
for (i in seq_len(nf)) {
  val <- unlist(lapply(get_res(files_nb[[i]]), function(x) mean(diag(x))))
  lines(sqd, val, pch = 19, col = pal[i], lwd = .7)
  points(sqd, val, pch = 19, col = pal[i], cex = .8)
}
addaxesD()
addlet("b")

# ml
plot(range(sqd), c(0.33, 1), type = "n", xlab = "", ylab = "")
addlet("c")
legend("bottomright", legend = idf, col = pal, pch = 19, bty = "n")

dev.off()

}