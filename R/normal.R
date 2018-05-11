
#' Confident effect sizes from from normal or t distributions
#'
#' A general purpose confident effect size function for where a normal or t distribution can be assumed. Calculates confident effect sizes based on a mean and standard deviation (normal distribution) or mean and scale (t distribution).
#'
#' @param mean A vector of means.
#'
#' @param sd A single or vector of standard deviations (or if t distribution, scales).
#'
#' @param df A single or vector of degrees of freedom, for t-distribution. Inf for normal distribution.
#'
#' @param fdr False Discovery Rate to control for.
#'
#' @param step Granularity of effect sizes to test.
#' 
#' @return
#'
#' See \code{\link{nest_confects}} for details of how to interpret the result.
#'
#'@export
normal_confects <- function(mean, sd=1, df=Inf, fdr=0.05, step=0.01) {
    n <- max(length(mean), length(sd), length(df))
    mean <- broadcast(mean, n)
    sd <- broadcast(sd, n)
    df <- broadcast(df, n)

    abs_mean <- abs(mean)

    pfunc <- function(indices, effect_size) {
        p_one <- function(i)
            sign_deciding_p(
                 abs_mean[i]/sd[i],
                 (abs_mean[i]-effect_size)/sd[i],
                 df[i])

        map_dbl(indices, p_one)
    }

    confects <- nest_confects(n, pfunc, fdr=fdr, step=step)
    effect <- mean[confects$table$index]
    confects$table$confect <- sign(effect) * confects$table$confect
    confects$table$effect <- effect

    confects
}