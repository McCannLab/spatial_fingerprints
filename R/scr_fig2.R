#' Figure 2
#'
#' Code to reproduce figure 2 (increase bio-tracers).
#'
#' @param path Path to results files for LDA and NBC.
#' @param file_ml Path to results file (a rds file) for MLP.
#'
#' @export

scr_fig2 <- function(path = "output/res_lda_nb/nbio/non_pca",
  file_ml = "output/res_f/res_combn_ml_nbio.rds") {

## helper functions
get_res <- function(file) {
  raw <- readRDS(file)
  lapply(raw, function(x) apply(x$mean, c(1, 2), mean))
}

get_res_alt <- function(file) {
  raw <- readRDS(file)
  lapply(raw, function(x) apply(x, c(1, 2), mean))
}

addlet <- function(let, x = 1) mtext(let, 3, at = x, cex = 1, font = 2)

add_ticks <- function() {
  axis(1, at = setdiff(sq_bt, c(5, 10, 15)), lwd = 0, lwd.ticks = .5,
  tck = -.025, labels = NA)
}

##
sqs <- c(1:3, 5, 10)
sq_bt <- seq_len(17)
msgInfo("Reading files for figure 2")
files_lda <- sprintf(paste0(path, "/nbio_lda_%02d.rds"), sqs)
files_nb <- sprintf(paste0(path, "/nbio_lda_%02d.rds"), sqs)
#
tmp <- readRDS(file_ml)
ml_reg <- lapply(tmp, function(x) aggregate(cbind(can, ru, us) ~ nbio, mean, data = x))


output_dir("output/figs")
msgInfo("Creating figure 2")

png("output/figs/fig2.png", width = 183, height = 121, units = "mm", res = 600)

par(mfrow = c(2, 3), las = 1, mar = c(2.5, 3.2, 2, .4), mgp = c(2.25, .6, 0))

## ROW 1 - performances
### P1 (LDA)
plot(range(sq_bt), c(0.33, 1), type = "n", xlab = "", ylab = "Performance")
pal <- gpuPalette("insileco")[c(2, 4, 3)]
for (i in 1:3) {
  val <- unlist(lapply(get_res(files_lda[[1]]), function(x) x[i, i]))
  lines(sq_bt, val, pch = 19, col = pal[i], lwd = .7)
  points(sq_bt, val, pch = 19, col = pal[i], cex = 1)
}
title(main = "LDA")
addlet("a")
add_ticks()
### P2 (NBC)
plot(range(sq_bt), c(0.33, 1), type = "n", xlab = "", ylab = "")
pal <- gpuPalette("insileco")[c(2, 4, 3)]
for (i in 1:3) {
  val <- unlist(lapply(get_res(files_nb[[1]]), function(x) x[i, i]))
  lines(sq_bt, val, pch = 19, col = pal[i], lwd = .7)
  points(sq_bt, val, pch = 19, col = pal[i], cex = 1)
}
title(main = "NBC")
addlet("b")
add_ticks()
### P3 (MLP)
plot(c(1,17), c(0.33, 1), type = "n", xlab = "", ylab = "")
for (i in 1:3) {
  val <- ml_reg[[1]][, i + 1]
  lines(sq_bt, val, pch = 19, col = pal[i], lwd = .7)
  points(sq_bt, val, pch = 19, col = pal[i], cex = 1)
}
legend("bottomright", legend = c("Canada", "Russia", "USA"), col = pal, pch = 19, bty = "n")
title(main = "MLP")
addlet("c")
add_ticks()

## ROW 2 - overall performances
par(mar = c(4, 3.2, .5, .4))

## P4 (LDA)
plot(c(1,17), c(0.33, 1), type = "n", xlab = "", ylab = "Overall performance")
pal <- colorRampPalette(c('grey10', 'grey70'))(5)
for (i in seq_len(5)) {
  val <- unlist(lapply(get_res(files_lda[[i]]), function(x) mean(diag(x))))
  lines(sq_bt, val, pch = 19, col = pal[i], lwd = .7)
  points(sq_bt, val, pch = 19, col = pal[i], cex = 1)
}
addlet("d")
add_ticks()
## P5 (NBC)
plot(c(1,17), c(0.33, 1), type = "n", xlab = "Number of bio-tracers combined", ylab = "")
pal <- colorRampPalette(c('grey10', 'grey70'))(5)
for (i in seq_len(5)) {
  val <- unlist(lapply(get_res(files_nb[[i]]), function(x) mean(diag(x))))
  lines(sq_bt, val, pch = 19, col = pal[i], lwd = .7)
  points(sq_bt, val, pch = 19, col = pal[i], cex = 1)
}
addlet("e")
add_ticks()
## P6 (MLP)
plot(range(sq_bt), c(0.33, 1), type = "n", xlab = "", ylab = "")
for (i in seq_len(5)) {
  id <- sqs[i]
  val <- apply(ml_reg[[id]][, 2:4], 1, mean)
  lines(sq_bt, val, pch = 19, col = pal[i], lwd = .7)
  points(sq_bt, val, pch = 19, col = pal[i], cex = 1)
}
addlet("f")
add_ticks()
legend("bottomright", legend = sqs, col = pal, pch = 19, bty = "n")

dev.off()


msgSuccess_fig("2", "output/figs")
invisible(0)

}