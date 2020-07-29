scr_fig3 <- function() {

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

addlet <- function(let) text(1, .98, let, cex = 1.2, font = 2)

addaxes <- function() {
  axis(1, at = seq(2, 16, 2), labels = NA, lwd = 0, lwd.ticks = .25, tck = -.025)
  axis(1, at = seq(1, 17, 2), labels = seq(1, 17, 2), lwd = 0, lwd.ticks = .5)
  axis(2,  lwd = 0,  lwd.ticks = .5)
  box()
}


files_lda <- sprintf("output/res_server_nb_lda/nbio_lda_%02d.rds",
  c(1))
files_nb <- sprintf("output/res_server_nb_lda/nbio_nb_%02d.rds",
  c(1))
# lda / pca / sample = 1
pca_lda <- unlist(lapply(get_res_pca("output/res_f/nbio_pca_lda_01.rds"), function(x) mean(diag(x))))
pca_nb <- unlist(lapply(get_res_pca("output/res_f/nbio_pca_nb_01.rds"), function(x) mean(diag(x))))


tmp <- readRDS('output/res_f/ml_nbio_pca.rds')
res_ml <- tmp[tmp$id_reg_test == tmp$id_reg_true, ]
ml_pca <- aggregate(prob~nbio, mean, data = res_ml)
tmp2 <- readRDS('output/res_f/res_ml_nbio.rds')
res_ml2 <- tmp2[tmp2$id_reg_test == tmp2$id_reg_true, ]
ml_reg <- aggregate(prob~nbio*id_comb, mean, data = res_ml2)
pca_ml <- split(ml_reg$prob, ml_reg$nbio)


png('output/figs/fig3.png',  width = 183, height = 70, units = "mm", res = 600)

par(mfrow = c(1, 3), las = 1, mar = c(4, 3.2, 1, .4), mgp = c(2.25, .6, 0))
plot(c(1, 17), c(.33, 1), type = "n", xlab = "", ylab = "overall performance", axes = FALSE)
boxplot(get_res_bb(files_lda[[1]]), col = "grey95", add = TRUE, pch = 19, border = "grey55", lwd = .8, cex = .5, axes = FALSE)
points(1:17, pca_lda, col = 1, pch = 19, cex = 1)
lines(1:17, pca_lda, col = 1, lwd = .7)
addaxes()


plot(c(1, 17), c(.33, 1), type = "n", xlab = "Number of biotracers", ylab = "", axes = FALSE)
boxplot(get_res_bb(files_nb[[1]]), col = "grey95", add = TRUE, pch = 19, border = "grey55", lwd = .8, cex = .5, axes = FALSE)
points(1:17, pca_nb, col = 1, pch = 19, cex = 1)
lines(1:17, pca_nb, col = 1, lwd = .7)
addaxes()

plot(c(1, 17), c(.33, 1), type = "n",  xlab = "", ylab = "", axes = FALSE)
boxplot(pca_ml, col = "grey95", add = TRUE, pch = 19, border = "grey55", lwd = .8, cex = .5, axes = FALSE)
points(1:17, ml_pca[, 2], col = 1, pch = 19, cex = 1)
lines(1:17, ml_pca[, 2], col = 1, lwd = .7)
addaxes()


dev.off()

}

