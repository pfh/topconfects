

#' General purpose function to find sets of discoveries for a range of effect sizes, controlling FDR
#'
#' Find sets of discoveries for a range of effect sizes, controlling the False Discovery Rate (FDR) for each set.
#'
#' @param n Number of items being tested.
#'
#' @param pfunc A function(indices, effect_size) to calculate p-values. Indices is a subset of 1:n giving the p-values to be computed.
#'
#' @param fdr False Discovery Rate to control for.
#'
#' @param step Granularity of effect sizes to test.
#'
#' @param full If TRUE, also include FDR-adjusted p-value that effect size is non-zero. Note that this is against the spirit of the topconfects approach.
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
#' The usage is as follows: To find a set of tests which have effect size greater than x with the specified FDR, take the rows with \code{abs(confect) >= x}. Once the set is selected, the confect values provide confidence bounds on the effect size with False Coverage-statement Rate (FCR) at the same level as the FDR.
#'
#' One may essentially take the top however many rows of the data frame and these will be the best set of results of that size to dependably have an effect size that is as large as possible. However if some genes have the same \code{abs(confect)}, all or none should be selected.
#'
#' Some rows in the output may be given the same \code{confect}, even if \code{step} is made small. This is an expected behaviour of the algorithm. (This is similar to FDR adjustment of p-values sometimes resulting in a run of the same adjusted p-value, even if all the input p-values are distinct.) 
#'
#' Some wrappers around this function may add a sign to the \code{confect} column, if it makes sense to do so. They will also generally add an \code{effect} column, containing an estimate of the effect size that aims to be unbiassed rather than a conservative lower bound.
#'
#' @examples
#'
#' # Find largest positive z-scores in a collection, 
#' # and place confidence bounds on them that maintain FDR 0.05.
#' z <- c(1,2,3,4,5)
#' pfunc <- function(i, effect_size) pnorm(z[i], mean=effect_size, lower.tail=FALSE)
#' nest_confects(length(z), pfunc, fdr=0.05)
#'
#' @export
nest_confects <- function(n, pfunc, fdr=0.05, step=0.001, full=FALSE) {
    indices <- seq_len(n)
    mags <- rep(NA, n)

    steps <- 0
    n_active <- n
    while(n_active > 0) {
        mag <- steps * step
        seq_active <- seq_len(n_active)
        p <- pfunc(indices[seq_active], mag)
        ordering <- order(p)
        p <- p[ordering]
        indices[seq_active] <- indices[ordering]
        mags[seq_active] <- mags[ordering]

        while(n_active > 0 && p[n_active] > fdr*n_active/n)
            n_active <- n_active - 1

        mags[seq_len(n_active)] <- mag
        steps <- steps + 1
    }

    table <- data.frame(rank=seq_len(n), index=indices, confect=mags)
    if (full) {
        table$fdr_zero <- 
            p.adjust(pfunc(seq_len(n), 0.0), method="BH")[ table$index ]
    }

    new("Topconfects", list(
        table=table,
        effect_desc = "effect size",
        fdr=fdr, step=step, pfunc=pfunc))
}


