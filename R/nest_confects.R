

#' General purpose function to find sets of discoveries for a range of effect sizes, controlling FDR
#'
#' Find sets of discoveries for a range of effect sizes, controlling the False Discovery Rate for each set.
#'
#' @param n Number of items being tested.
#'
#' @param pfunc A function(indices, effect_size) to calculate p-values. Indices is a subset of 1:n giving the p-values to be computed.
#'
#' @param fdr False Discovery Rate to control for.
#'
#' @param step Granularity of effect sizes to test.
#'
#' @return A "Topconfects" object, containing a table of results and various associated information.
#'
#' The most important part of this object is the $table element, a data frame with the following columns:
#'
#' \itemize{
#'     \item \code{rank} - Ranking by \code{confect} and for equal \code{confect} by p-value at that effect size.
#'     \item \code{index} - Number of the test, between 1 and n.
#'     \item \code{confect} - CONfident efFECT size.
#' }
#'
#' The usage is as follows: To find a set of tests which have effect size at least \code{x} with the specified FDR, take the rows with \code{abs(confect) >= x}.
#'
#' Some tests may have been given the same \code{confect}. To maintain the FDR, all or none should be chosen.
#'
#' With this caveat understood, one may essentially take the top however many rows of the data frame and these will be the best set of results of that size to dependably have an effect size that is as large as possible.
#'
#' Some wrappers around this function may add a sign to the \code{confect} column, if it makes sense to do so. They will also generally add an \code{effect} column, containing an estimate of the effect size that aims to be unbiassed rather than a conservative lower bound.
#'
#' @export
nest_confects <- function(n, pfunc, fdr=0.05, step=0.01) {
    indices <- seq_len(n)
    mags <- rep(NA, n)

    mag <- 0.0
    n_active <- n
    while(n_active > 0) {
        seq_active <- seq_len(n_active)
        p <- pfunc(indices[seq_active], mag)
        ordering <- order(p)
        p <- p[ordering]
        indices[seq_active] <- indices[ordering]
        mags[seq_active] <- mags[ordering]

        while(n_active > 0 && p[n_active] > fdr*n_active/n)
            n_active <- n_active - 1

        mags[seq_len(n_active)] <- mag
        mag <- mag + step
    }

    new("Topconfects", list(
        table=data.frame(rank=seq_len(n), index=indices, confect=mags),
        effect_desc = "effect size",
        fdr=fdr, step=step, pfunc=pfunc))
}


