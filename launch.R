devtools::load_all()
i <- as.numeric(commandArgs(trailingOnly = TRUE))
mx_cb <- 500
n_rep <- 200
n_rep2 <- 10000

output_dir("res/")
mysave <- function(obj, scr, meth, i, path = "res/") {
  fl <-  sprintf("%s%s_%s_%02d.rds", path, scr, meth, i)
  saveRDS(obj, file = fl)
  invisible(fl)
}

if (i < 11) {
  k <- i
  out <- simu_nbio("lda", mxcb = mx_cb, nrep = n_rep, nsample = k, ndistr = 20)
  mysave(out, "nbio", "lda", k)
} else if (i > 10 & i < 21) {
  k <- i - 10
  out <- simu_nbio("nb", mxcb = mx_cb, nrep = n_rep, nsample = k, ndistr = 20)
  mysave(out, "nbio", "nb", k)
} else if (i > 20 & i < 31) {
  k <- i - 20
  out <- simu_nbio_order(method = "lda", nrep = n_rep2, nsample = k, ndistr = 20)
  mysave(out, "nbio_pca", "lda", k)
} else if (i > 30 & i < 41) {
  k <- i - 30
  out <- simu_nbio_order(method = "nb", nrep = n_rep2, nsample = k, ndistr = 20)
  mysave(out, "nbio_pca", "nb", k)
} else if (i > 40 & i < 51) {
  k <- i - 40
  out <- simu_noise("lda", mxcb = mx_cb, nrep = n_rep, nsample = 1, nbio = k, ndistr = 20)
  mysave(out, "noise", "lda", k)
} else if (i > 50 & i < 61) {
  k <- i - 50
  out <- simu_noise("nb", mxcb = mx_cb, nrep = n_rep, nsample = 1, nbio = k, ndistr = 20)
  mysave(out, "noise", "nb", k)
} else if (i > 60 & i < 71) {
  k <- i - 60
  out <- simu_nsample("lda", mxcb = mx_cb, nrep = n_rep, ndistr = 15, nbio = k)
  mysave(out, "nsample", "lda", k)
} else if (i > 70 & i < 81) {
  k <- i - 70
  out <- simu_nsample("nb", mxcb = mx_cb, nrep = n_rep, ndistr = 15, nbio = k)
  mysave(out, "nsample", "nb", k)
} else if (i > 80 & i < 91) {
  k <- i - 80
  out <- simu_ndistr("lda", mxcb = mx_cb, nrep = n_rep, nsample = 1, nbio = k)
  mysave(out, "ndistr", "lda", k)
} else if (i > 90 & i < 101) {
  k <- i - 90
  out <- simu_ndistr("nb", mxcb = mx_cb, nrep = n_rep, nsample = 1, nbio = k)
  mysave(out, "ndistr", "nb", k)
}

