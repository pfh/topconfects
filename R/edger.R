
#' Confident effect sizes based on the edgeR Quasi-Likelihood method, both linear and non-linear
#'
#' For all possible absolute log2 fold changes, which genes have at least this fold change at a specified False Discovery Rate?
#'
#' @param fit An edgeR DGEGLM object produced using \code{glmQLFit}.
#'
#' @param coef Coefficient to test, as per \code{glmTreat}. Use either coef or contrast or effect.
#'
#' @param contrast Contrast to test, as per \code{glmTreat}. Use either coef or contrast or effect.
#'
#' @param effect A non-linear effect, created with one of the effect_... functions. Use either coef or contrast or effect.
#'
#' @param fdr False Discovery Rate to control for.
#'
#' @param step Granularity of log2 fold changes to test.
#'
#' @param null "null" parameter passed through to edger::glmTreat (if coef or contrast given). Choices are "worst.case" or "interval". Note that the default here is "worst.case", to be consistent with other functions in topconfects. This differs from the default for glmTreat.
#'
#' @return
#'
#' See \code{\link{nest_confects}} for details of how to interpret the result.
#'
#' Technical note: when using a non-linear effect size: Signed confects are based on TREAT-style p-values. Unsigned confects (generally with df>1) are based on comparing the best fit within the H0 region to the best fit overall, which may up to double p-values.
#'
#' @export
edger_confects <- function(fit, coef=NULL, contrast=NULL, effect=NULL, fdr=0.05, step=0.01, null="worst.case") {
    assert_that(is(fit, "DGEGLM"))
    assert_that((!is.null(coef)) + (!is.null(contrast)) + (!is.null(effect)) == 1)

    if (!is.null(effect))
        confects <- edger_nonlinear_confects(fit, effect, fdr=fdr, step=step)
    else {
        n <- nrow(fit)

        top_tags <-
            glmQLFTest(fit,coef=coef,contrast=contrast) %>%
            topTags(n=n,sort.by="none")

        pfunc <- function(i, mag) {
            if (mag == 0.0)
                top_treats <- top_tags
            else
                top_treats <-
                    glmTreat(fit, coef=coef, contrast=contrast, lfc=mag, null=null) %>%
                    topTags(n=n, sort.by="none")

            top_treats$table$PValue[i]
        }

        confects <- nest_confects(n, pfunc, fdr=fdr, step=step)
        confects$effect_desc <- "log2 fold change"
        logFC <- top_tags$table$logFC[confects$table$index]
        confects$table$confect <- sign(logFC) * confects$table$confect
        confects$table$effect <- logFC
        confects$table$logCPM <- top_tags$table$logCPM[confects$table$index]
        confects$table$name <- rownames(top_tags$table)[confects$table$index]
    }

    confects$edger_fit <- fit

    if (!is.null(fit$genes)) {
        confects$table <- cbind(confects$table, fit$genes[confects$table$index,])
    }

    confects
}
