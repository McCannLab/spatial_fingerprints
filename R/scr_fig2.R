scr_fig2 <- function() {

get_res <- function(file) {
  raw <- readRDS(file)
  lapply(raw, function(x) apply(x$mean, c(1, 2), mean))
}

get_res_alt <- function(file) {
  raw <- readRDS(file)
  lapply(raw, function(x) apply(x, c(1, 2), mean))
}

addlet <- function(let, x = 1, y = .98) text(x, y, let, cex = 1.2, font =2)


sqs <- c(1:3, 5, 10)
files_lda <- sprintf("output/res_server_nb_lda/nbio_lda_%02d.rds", sqs)
files_nb <- sprintf("output/res_server_nb_lda/nbio_nb_%02d.rds", sqs)
#
tmp <- readRDS('output/res_f/res_combn_ml_nbio.rds')
ml_reg <- lapply(tmp, function(x) aggregate(cbind(can, ru, us) ~ nbio, mean, data = x))



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
  val <- ml_reg[[1]][, i + 1]
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
addlet("d", y = .34)

# nb
plot(c(1,17), c(0.33, 1), type = "n", xlab = "Number of biotracers", ylab = "")
pal <- colorRampPalette(c('grey10', 'grey70'))(5)
for (i in 1:5) {
  val <- unlist(lapply(get_res(files_nb[[i]]), function(x) mean(diag(x))))
  lines(1:17, val, pch = 19, col = pal[i], lwd = .7)
  points(1:17, val, pch = 19, col = pal[i], cex = 1)
}
addlet("e", y = .34)

# ml
plot(c(1,17), c(0.33, 1), type = "n", xlab = "", ylab = "")

for (i in 1:5) {
  id <- c(1:3, 5, 10)[i]
  val <- apply(ml_reg[[id]][, 2:4], 1, mean)
  lines(1:17, val, pch = 19, col = pal[i], lwd = .7)
  points(1:17, val, pch = 19, col = pal[i], cex = 1)
}
addlet("f", y = .34)
legend("bottomright", legend = sqs, col = pal, pch = 19, bty = "n")

dev.off()

}