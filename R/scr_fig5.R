scr_fig5 <- function(file = "output/res_lda_nb/all123/res_lda_123.rds") {


tmp <- lapply(readRDS(file),
  function(x) apply(x$mean, 3, function(y) mean(diag(y))))

hh <- get_data_ready()
res_ind <- data.frame(
  perf_ind = tmp[[1]],
  var = 0,
  var_wtn = 0,
  var_btw = 0,
  mean_distb = 0  # distance between centroids
)

add_ticks <- function(at) axis(1, at = at, lwd = 0, lwd.ticks = .5, tck = -.025, labels = NA)



# intertian inter / intra + distance between points
for (i in 3:19) {
  res_ind$var[i - 2] <- sum(hh[, i]^2)
  val1 <- hh[1:30, i]
  val2 <- hh[31:60, i]
  val3 <- hh[61:90, i]
  res_ind$var_wtn[i - 2] <- myvar(val1) + myvar(val2) + myvar(val3)
  m1 <- mean(hh[1:30, i])
  m2 <- mean(hh[31:60, i])
  m3 <- mean(hh[61:90, i])
  res_ind$var_btw[i - 2] <- 30*(m1^2 + m2^2 + m3^2)
  res_ind$mean_distb[i - 2] <- mean((m1-m2)^2, (m1-m3)^2, (m2-m3)^2)
}

# Same work for pair of bio-tracers
comb_bio2 <- data.frame(t(combn(3:19, 2)), perf = 0)
comb_bio2$max_perf1 <- apply(comb_bio2, 1, function(x) getmax(x[1:2], perf = res_ind$perf_ind))
comb_bio2$mean_perf2 <- apply(comb_bio2, 1, function(x) getsum(x[1:2], perf = res_ind$perf_ind))
comb_bio2$categ <- apply(comb_bio2, 1, function(x) get_categ2(x[1], x[2]))
comb_bio2$mean_distb <- apply(comb_bio2, 1, function(x) getdist(hh[x[1:2]]))
comb_bio2 <- cbind(comb_bio2, t(apply(comb_bio2, 1, function(x) inertia(hh[x[1:2]]))))

names(comb_bio2)[1:2] <- paste0("bio", 1:2)
comb_bio2$perf <- tmp[[2]]


getval <- function(ids3, FUN = max) {
  FUN(
    comb_bio2$perf[
    apply(comb_bio2[,1:2], 1, function(x) sum(x %in% ids3) == 2)
    ]
  )
}

comb_bio3 <- data.frame(t(combn(3:19, 3)), perf = 0)
names(comb_bio3)[1:3] <- paste0("bio", 1:3)
comb_bio3$max_perf2 <- apply(comb_bio3, 1, function(x) getval(x[1:3]))
comb_bio3$mean_perf2 <- apply(comb_bio3, 1, function(x) getval(x[1:3], FUN = mean))
comb_bio3$mean_distb <- apply(comb_bio3, 1, function(x) getdist(hh[x[1:3]]))
comb_bio3 <- cbind(comb_bio3, t(apply(comb_bio3, 1, function(x) inertia(hh[x[1:3]]))))

comb_bio3$perf <- tmp[[3]]


colr <- "#f63267"

png("output/figs/fig5f.png", width = 130, height = 60, units = "mm", res = 600)
layout(matrix(1:3, ncol = 3), widths = c(1, 1, .35))
par(las = 1, mar = c(3.2, 3.8, 1, 0), mgp = c(2.6, .7, 0), cex.axis = .8, cex.lab = .9)

vc_pos <- rep(3, 17)
vc_pos[c(8, 13, 9)] <- 1
vc_pos[c(3, 14, 16)] <- 2
vc_pos[c(6, 11, 15, 13)] <- 4
rgy <- c(0.3, 0.6)
 # plot(apply(abs(res_ind[,-1]), 1, mean), res_ind[,1], pch = 19)
palg <- c(rep("grey75", 3), rep("grey25", 14))
dpalg <- sapply(palg, darken, 25)
#
nmb <- gsub("n", "n-", gsub("_", ":", gsub("FA_", "", names(hh))))
#
plot(100*res_ind$var_btw/89, res_ind$perf_ind, bg = palg, col = dpalg, pch = 21,
  cex = .9, ylim = rgy, xlim = 100*c(0, .55), xlab = "", ylab = "Performance")
par(mgp = c(2, .7, 0))
title(xlab = "Inter-regions inertia (%)")
text(100*res_ind$var_btw/89, res_ind$perf_ind, 1:17, pos = vc_pos, offset = .3, col = palg, cex = .8)
add_ticks(seq(5, 55, 5))
mtext("a", 3, at = 0, font = 2, cex = .8)
##
vc_pos[c(3, 6, 10)] <- 1
vc_pos[c(1)] <- 4
vc_pos[c(8:9, 16)] <- 3
plot(res_ind$mean_distb, res_ind$perf_ind, bg = palg, col = dpalg, pch = 21,
  ylim = rgy, xlab = "", ylab = "", cex = .9)
title(xlab = "Mean distance between centroids")
text(res_ind$mean_distb, res_ind$perf_ind, 1:17, pos = vc_pos, offset = .3, col = palg, cex = .8)
add_ticks(seq(.25, 3.25, .5))
mtext("b", 3, at = 0, font = 2, cex = .8)


par(mar = c(4.6, .4, 1, .8))
plot(c(0, 1), c(1, 17), type = "n", axes = FALSE, ann = FALSE)
for (i in seq_len(17)) {
  text(.44, 18 - i, paste0(i, ". "), col = dpalg[i], pos = 2, cex = .8)
  text(.2, 18 - i, nmb[i + 2], col = dpalg[i], pos = 4, cex = .8)
}

dev.off()



png("output/figs/fig6.png", width = 89, height = 52, units = "mm", res = 600)


par(mfrow = c(1, 2), yaxs = "i", las = 1, mar = c(4, 4, 1.5, .5), mgp = c(2.5, .7, 0), cex = .45)
plot(comb_bio2$max_perf1, comb_bio2$perf, pch = 20, col = "grey10",
  xlab = c("Best single biotracer performance"),
  ylab = c("Performance of a pair of bio-tracers"), xlim = c(.3, .72), ylim = c(.3, .85), cex = .6)
abline(a = 0, b = 1, lty = 3, col = colr, lwd = 1)
mtext("a", 3, at = 0.3, font = 2, cex = .7)

plot(comb_bio3$max_perf2, comb_bio3$perf, pch = 20, col = "grey10",
  xlab = "Performance of the best pair of bio-tracers",
  ylab = "Performance of a combination (n=3)",
  xlim = c(.3, .72) , ylim = c(.3, .85), cex = .6)
abline(a = 0, b = 1, lty = 3, col = colr, lwd = 1)
mtext("b", 3, at = 0.3, font = 2, cex = .7)

dev.off()





png("output/figs/fig_7b.png",  width = 89, height = 52, units = "mm", res = 600)

add_rsq <- function(rsq, x =7.2, y= .32) {
  text(x, y, cex = .9, pos = 2, labels = expression(R^2==""), offset = 0)
  text(x, y, cex = .9, pos = 4,
      labels = paste0(format(100*rsq, digit = 3), "%"), offset = .1)
}

par(mfrow = c(1, 2), yaxs = "i", las = 1, mar = c(4, 4, 1.5, .5), mgp = c(2.6, .7, 0), cex = .45)

plot(comb_bio2$mean_distb, comb_bio2$perf, pch = 20, col = "grey10",
  xlab = "Mean distance between centroids",
  ylab = "Performance of a pair of biotracer", ylim = c(.3, .82), cex = .9)
f <- fitexp(comb_bio2, lty = 2, col = colr,  lwd = 1)
rsq <- 1 - deviance(f)/deviance(lm(perf~1, data = comb_bio2))
add_rsq(rsq, x = 5)
mtext("a", 3, at = 0, font = 2, cex = .7)

plot(comb_bio3$mean_distb, comb_bio3$perf, pch = 20, col = "grey10",
  xlab = "Mean distance between centroids",
  ylab = "Performance of a combination (n=3)", ylim = c(.3, .82), cex = .9)
f <- fitexp(comb_bio3, lty = 2, col = colr,  lwd = 1)
rsq <- 1 - deviance(f)/deviance(lm(perf~1, data = comb_bio3))
add_rsq(rsq)
mtext("b", 3, at = 0, font = 2, cex = .7)


dev.off()


png("output/figs/fig7.png", width = 89, height = 50, units = "mm", res = 600)

par(mfrow = c(1, 2), yaxs = "i", las = 1, mar = c(4, 4, 1.5, .5), mgp = c(2.6, .65, 0), cex = .45)

prop2 <- comb_bio2$inert_btw/comb_bio2$inert_tot
plot(100*prop2, comb_bio2$perf, pch = 20, col = "grey10",
  xlab = "",
  ylab = "Performance of a pair of biotracer", ylim = c(.3, .82), cex = .8)
f <- fitexp2(comb_bio2, 100*prop2, lty = 2, col = colr,  lwd = 1)
rsq <- 1 - deviance(f)/deviance(lm(perf~1, data = comb_bio2))
add_rsq(rsq, x = 45)
title(xlab = "Inter-regions inertia (%)")
mtext("a", 3, at = 0, font = 2, cex = .7)
add_ticks(seq(5, 50, 5))

par(mgp = c(2.6, .65, 0))
prop3 <- comb_bio3$inert_btw/comb_bio3$inert_tot
plot(100*prop3, comb_bio3$perf, pch = 20, col = "grey10",
  xlab = "",
  ylab = "Performance of a combination (n=3)", ylim = c(.3, .82), cex = .8)
f <- fitexp2(comb_bio3, 100*prop3, lty = 2, col = colr,  lwd = 1)
rsq <- 1 - deviance(f)/deviance(lm(perf~1, data = comb_bio3))
add_rsq(rsq, x = 45)
title(xlab = "Inter-regions inertia (%)")
mtext("b", 3, at = 0, font = 2, cex = .7)
add_ticks(seq(5, 50, 5))

dev.off()


}
