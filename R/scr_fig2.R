scr_fig2 <- function() {

get_res <- function(file) {
  raw <- readRDS(file)
  lapply(raw, function(x) apply(x$mean, c(1, 2), mean))
}

get_res_sd <- function(file) {
  raw <- readRDS(file)
  lapply(raw, function(x) apply(x$sd, c(1, 2), mean))
}


get_res_alt <- function(file) {
  raw <- readRDS(file)
  lapply(raw, function(x) apply(x, c(1, 2), mean))
}

addlet <- function(let) text(1, .98, let, cex = 1.2, font =2)


files_lda <- sprintf("output/res_server_nb_lda/nbio_lda_%02d.rds",
  c(1:3, 5, 10))
files_nb <- sprintf("output/res_server_nb_lda/nbio_nb_%02d.rds",
  c(1:3, 5, 10))
#
tmp <- readRDS('output/res_f/res_ml_nbio.rds')
res_ml <- tmp[tmp$id_reg_test == tmp$id_reg_true, ]
ml_reg <- aggregate(prob~nbio*id_reg_true, mean, data = res_ml)
ml_sam <- aggregate(prob~nbio, mean, data = res_ml)



png("output/figs/fig2.png", width = 183, height = 121, units = "mm", res = 600)

# FIGURE 1
par(mfrow = c(2, 3), las = 1, mar = c(2, 3.2, 2, .4), mgp = c(2.25, .6, 0))

## ROW 1

plot(c(1, 17), c(0.33, 1), type = "n", xlab = "", ylab = "Performance")
pal <- gpuPalette("insileco")[c(2, 4, 3)]
for (i in 1:3) {
  val <- unlist(lapply(get_res(files_lda[[1]]), function(x) x[i, i]))
  lines(1:17, val, pch = 19, col = pal[i], lwd = .7)
  points(1:17, val, pch = 19, col = pal[i], cex = 1)
}
title(main = "LDA")
addlet("a")

plot(c(1, 17), c(0.33, 1), type = "n", xlab = "", ylab = "")
pal <- gpuPalette("insileco")[c(2, 4, 3)]
for (i in 1:3) {
  val <- unlist(lapply(get_res(files_nb[[1]]), function(x) x[i, i]))
  lines(1:17, val, pch = 19, col = pal[i], lwd = .7)
  points(1:17, val, pch = 19, col = pal[i], cex = 1)
}
title(main = "NB")
addlet("b")

# ml
plot(c(1,17), c(0.33, 1), type = "n", xlab = "", ylab = "")
for (i in 1:3) {
  val <- ml_reg$prob[ml_reg$id_reg == i]
  lines(1:17, val, pch = 19, col = pal[i], lwd = .7)
  points(1:17, val, pch = 19, col = pal[i], cex = 1)
}
legend("bottomright", legend = c("Canada", "Russia", "USA"), col = pal, pch = 19, bty = "n")
title(main = "ML")
addlet("c")


## ROW 2
par(mar = c(4, 3.2, 0, .4))

## lda
plot(c(1,17), c(0.33, 1), type = "n", xlab = "", ylab = "Overall performance")
pal <- colorRampPalette(c('grey10', 'grey70'))(5)
for (i in 1:5) {
  val <- unlist(lapply(get_res(files_lda[[i]]), function(x) mean(diag(x))))
  lines(1:17, val, pch = 19, col = pal[i], lwd = .7)
  points(1:17, val, pch = 19, col = pal[i], cex = 1)
}
addlet("d")

# nb
plot(c(1,17), c(0.33, 1), type = "n", xlab = "Number of biotracers", ylab = "")
pal <- colorRampPalette(c('grey10', 'grey70'))(5)
for (i in 1:5) {
  val <- unlist(lapply(get_res(files_nb[[i]]), function(x) mean(diag(x))))
  lines(1:17, val, pch = 19, col = pal[i], lwd = .7)
  points(1:17, val, pch = 19, col = pal[i], cex = 1)
}
addlet("e")

# ml
plot(c(1,17), c(0.33, 1), type = "n", xlab = "", ylab = "")
points(ml_sam[, 1], ml_sam[, 2], pch = 19)
lines(ml_sam[, 1], ml_sam[, 2], pch = 19)
addlet("f")

#
legend("bottomright", legend = 1:5, col = pal, pch = 19, bty = "n")

dev.off()

}