
#' Confident effect sizes from from normal or t distributions
#'
#' A general purpose confident effect size function for where a normal or t distribution of errors can be assumed. Calculates confident effect sizes based on an estimated effect and standard deviation (normal distribution), or mean and scale (t distribution).
#'
#' @param effect A vector of estimated effects.
#'
#' @param se A single or vector of standard deviations (or if t distribution, scales).
#'
#' @param df A single or vector of degrees of freedom, for t-distribution. Inf for normal distribution.
#'
#' @param fdr False Discovery Rate to control for.
#'
#' @param step Granularity of effect sizes to test.
#'
#' @param full Include some further statistics used to calculate confects in the output.
#' 
#' @return
#'
#' See \code{\link{nest_confects}} for details of how to interpret the result.
#'
#'@export
normal_confects <- function(effect, se, df=Inf, fdr=0.05, step=0.001, full=FALSE) {
    n <- max(length(mean), length(se), length(df))
    effect <- broadcast(effect, n)
    se <- broadcast(se, n)
    df <- broadcast(df, n)

    abs_effect <- abs(effect)

    pfunc <- function(indices, mag) {
        abs_effect_i <- abs_effect[indices]
        se_i <- se[indices]
        df_i <- df[indices]

        tstat_right <- (abs_effect_i-mag)/se_i
        tstat_left <- (abs_effect_i+mag)/se_i

        pt(tstat_right, df=df_i, lower.tail=FALSE) + 
            pt(tstat_left,df=df_i, lower.tail=FALSE)
    }

    confects <- nest_confects(n, pfunc, fdr=fdr, step=step)
    ranked_effect <- effect[confects$table$index]
    confects$table$confect <- sign(ranked_effect) * confects$table$confect
    confects$table$effect <- ranked_effect

    if (full) {
        confects$table$se <- se[confects$table$index]
        confects$table$df <- df[confects$table$index]
    }

    confects
}