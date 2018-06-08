
#' Confident effect sizes from from normal or t distributions
#'
#' A general purpose confident effect size function for where a normal or t distribution of errors can be assumed. Calculates confident effect sizes based on an estimated effect and standard deviation (normal distribution), or mean and scale (t distribution).
#'
#' @param effect A vector of estimated effects.
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
normal_confects <- function(effect, sd, df=Inf, fdr=0.05, step=0.01) {
    n <- max(length(mean), length(sd), length(df))
    effect <- broadcast(effect, n)
    sd <- broadcast(sd, n)
    df <- broadcast(df, n)

    abs_effect <- abs(effect)

    pfunc <- function(indices, mag) {
        abs_effect_i <- abs_effect[indices]
        sd_i <- sd[indices]
        df_i <- df[indices]

        tstat_right <- (abs_effect_i-mag)/sd_i
        tstat_left <- (abs_effect_i+mag)/sd_i

        pt(tstat_right, df=df_i, lower.tail=FALSE) + 
            pt(tstat_left,df=df_i, lower.tail=FALSE)
    }

    confects <- nest_confects(n, pfunc, fdr=fdr, step=step)
    ranked_effect <- effect[confects$table$index]
    confects$table$confect <- sign(ranked_effect) * confects$table$confect
    confects$table$effect <- ranked_effect

    confects
}