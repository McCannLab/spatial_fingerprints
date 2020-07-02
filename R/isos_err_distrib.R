#' Sample in a distribution of values for which an error value in associated.
#'
#' Sample in a distribution of values for which an error value in associated.
#'
#' @param vc_values List of values to be sampled.
#' @param vc_errors List of error values to be sampled.
#' @param sample_size size of the samples to be drawn.
#'
#' @importFrom stats rnorm
#' @export
#'
#' @return
#' A vector of size `sample_size`
#'
#' @examples
#' isos_err_distrib(1:10, rep(2,10), 100)

isos_err_distrib <- function(vc_values, vc_errors, sample_size) {
    stopifnot(length(vc_values) == length(vc_errors))
    sz <- length(vc_values)
    id <- sample(sz, sample_size, replace = TRUE)
    unlist(lapply(id, function(x) rnorm(1, vc_values[x], vc_errors[x])))
}
