#' @name isos_simulations
#' @title simulations for isoscape
#'
#' @description Functions to compute probability of origins to a given sample.
#'
#' @param ls_values List of values to be sampled.
#' @param sample_size size of the samples to be drawn.
#' @param sample values of a sample.
#' @param ls_density a list of object of class \code{\link[stats]{density}}.
#' @param dens1 an integer that identities the first distribution to be mixed (in \code{ls_density}).
#' @param dens2 an integer that identities the second distribution to be mixed (in \code{ls_density}).
#' @param perc a vector of numeric indicating the percentage of dens1 included in the sample (assuming the remaining is made of dens2).
#' @param density an object of class \code{\link[stats]{density}}.
#' @param nrep an integer indicating the number of repetition (see `details` section).
#'
#' @details
#' For every density included in `ls_density` the likelihood that the sample is
#' coming from a specific region is computed. Then likelihood values are
#' compared to determine the probability for `sample` of having been drawn
#' from the density (there is a reason why the ratio work). Probability for
#' all densities are returned and the corresponding vector is ordered following
#' `ls_density`.
#'
#' @export
#'
#' @return
#' A matrix

isos_simulations <- function(ls_values, sample_size, ls_density, nrep) {
    # check
    stopifnot(length(ls_density) == length(ls_values))
    #
    lg <- length(ls_density)
    out <- matrix(0, lg, lg)
    #
    for (k in seq_len(nrep)) {
        for (i in seq_len(lg)) {
            out[i, ] <- out[i, ] + isos_sample(sample(ls_values[[i]], sample_size,
                replace = TRUE), ls_density)
        }
    }
    # average over the number of repetitions
    out/nrep
}


#' @describeIn isos_simulations Analyse for one sample.
#' @return
#' A vector of probabilities.
isos_simulations_mixture2 <- function(ls_values, sample_size, ls_density, dens1,
    dens2, perc, nrep) {
    # check
    stopifnot(length(ls_density) == length(ls_values))
    stopifnot(perc <= 100 & perc >= 0)
    #
    lg <- length(perc)
    out <- matrix(0, lg, length(ls_density))
    #
    for (k in seq_len(nrep)) {
        for (i in seq_len(lg)) {
            lg2 <- length(ls_values[[dens1]])
            smpl <- sample(ls_values[[dens1]], sample_size, replace = TRUE)
            nbr <- floor(perc[i]/100 * lg2)
            if (nbr)
                smpl[sample(1:lg2, nbr)] <- sample(ls_values[[dens2]], nbr, replace = TRUE)
            out[i, ] <- out[i, ] + isos_sample(sample(smpl, sample_size, replace = TRUE),
                ls_density)
        }
    }
    # average over the number of repetitions
    out/nrep
}



#' @describeIn isos_simulations Analyse for one sample.
#' @return
#' A vector of probabilities of origin for every distribution in ls_density.
isos_sample <- function(sample, ls_density) {
    # check
    stopifnot(class(ls_density) == "list")
    #
    lg <- length(ls_density)
    tmp <- numeric(lg)
    for (i in seq_len(lg))
        tmp[i] <- isos_likelihood(sample, ls_density[[i]])
    get_proba(tmp)
}


#' @describeIn isos_simulations Compute the likelihood associated with one sample.
#' @return A likelihood value.
isos_likelihood <- function(sample, density) {
    # check
    stopifnot(class(density) == "density")
    #
    ls_tmp <- lapply(sample, FUN = get_loglik, density = density)
    sum(-log(unlist(ls_tmp)))
}


get_loglik <- function(x, density) {
    val <- (x - density$x)^2
    density$y[which.min(val)]
}

get_proba <- function(vec) {
    # check
    if (all(is.infinite(vec)))
        return(rep(1/length(vec), length(vec)))
    #
    vec <- vec - vec[which.min(vec)]
    # 0 and infinite are special cases handled separately
    id <- vec > 0 & !is.infinite(vec)
    #
    vec[vec == 0] <- 1
    vec[is.infinite(vec)] <- 0
    if (length(id))
        vec[id] <- exp(-vec[id])
    #
    vec/sum(vec)
}
