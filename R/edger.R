
#' Confident log2 fold change based on the edgeR Quasi-Likelihood method
#'
#' For all possible absolute log2 fold changes (LFC), which genes have at least this fold change at a specified False Discovery Rate?
#'
#' Results are presented in a table such that for any given LFC, if the reader chooses the genes with abs(confect) less than this they are assured that this set of genes has at least this LFC (with the specified FDR). The confect column may also be viewed as a confidence bound on the LFC of each gene, with a dynamic correction for multiple testing.
#'
#' @param fit An edgeR DGEGLM object produced using \code{glmQLFit}.
#'
#' @param coef Coefficient to test, as per \code{glmTreat}. Use either coef or contrast or effect.
#'
#' @param contrast Contrast to test, as per \code{glmTreat}. Use either coef or contrast or effect.
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
#' @export
edger_confects <- function(fit, coef=NULL, contrast=NULL, fdr=0.05, step=0.01, null="worst.case") {
    assert_that(is(fit, "DGEGLM"))
    assert_that((!is.null(coef)) + (!is.null(contrast)) == 1)

    n <- nrow(fit)

    tested <- edgeR::glmQLFTest(fit,coef=coef,contrast=contrast)
    top_tags <- edgeR::topTags(tested, n=n,sort.by="none")

    pfunc <- function(i, mag) {
        if (mag == 0.0)
            top_treats <- top_tags
        else {
            treat_tested <- edgeR::glmTreat(fit, coef=coef, contrast=contrast, lfc=mag, null=null)
            top_treats <- edgeR::topTags(treat_tested, n=n, sort.by="none")
        }

        top_treats$table$PValue[i]
    }

    confects <- nest_confects(n, pfunc, fdr=fdr, step=step)
    confects$effect_desc <- "log2 fold change"
    logFC <- top_tags$table$logFC[confects$table$index]
    confects$table$confect <- sign(logFC) * confects$table$confect
    confects$table$effect <- logFC
    confects$table$logCPM <- top_tags$table$logCPM[confects$table$index]
    confects$table$name <- rownames(top_tags$table)[confects$table$index]

    confects$edger_fit <- fit

    if (!is.null(fit$genes)) {
        confects$table <- cbind(confects$table, fit$genes[confects$table$index,,drop=FALSE])
    }

    confects
}
