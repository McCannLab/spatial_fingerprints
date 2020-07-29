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
  k = i
  out <- simu_ndistr("lda", mxcb = mx_cb, nrep = n_rep, nsample = 1, nbio = k, ndistr = 20)
  mysave(out, "noise", "lda", k)
} else if (i == 11) {
  k = 15
  out <- simu_ndistr("nb", mxcb = mx_cb, nrep = n_rep, nsample = 1, nbio = k, ndistr = 20)
  mysave(out, "noise", "nb", k)
} else if (i == 12) {
  k = 15
  out <- simu_ndistr("lda", mxcb = mx_cb, nrep = n_rep, nsample = 1, nbio = k)
  mysave(out, "ndistr", "lda", k)
}
