

#' Find sets of discoveries for a range of effect sizes, controlling the False Discovery Rate for each set.
#'
#' @param n Number of items being tested.
#'
#' @param pfunc A function(indices, effect_size) to calculate p-values. Indices is a subset of 1:n giving the p-values to be computed.
#'
#' @param fdr False Discovery Rate to control for.
#'
#' @param max Maximum effect size to test for.
#'
#' @param step Granularity of effect sizes to test.
#'
#' @export
nest_confects <- function(n, pfunc, fdr=0.01, max=30.0, step=0.05) {
    indices <- seq_len(n)
    mags <- rep(NA, n)

    mag <- 0.0
    n_active <- n
    while(n_active > 0 && mag <= max) {
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

    data.frame(rank=seq_len(n), index=indices, confect=mags)
}


