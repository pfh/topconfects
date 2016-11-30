
#' Confident log2 fold changes using edgeR's \code{glmTreat}
#'
#' For all possible absolute log2 fold changes, which genes have at least this fold change at a specified False Discovery Rate?
#'
#' @param fit An edgeR DGEGLM object.
#'
#' @param coef Coefficient to test, as per \code{glmTreat}. Use either coef or contrast.
#'
#' @param contrast Contrast to test, as per \code{glmTreat}. Use either coef or contrast.
#'
#' @param fdr False Discovery Rate to control for.
#'
#' @param max Maximum log2 fold change to test for.
#'
#' @param step Granularity of log2 fold changes to test.
#'
#' @return
#'
#' See \code{\link{nest_confects}} for details of how to interpret the result.
#'
#' @export
edger_confects <- function(fit, coef=ncol(fit$design), contrast=NULL, fdr=0.05, max=30.0, step=0.05) {
    assert_that(is(fit, "DGEGLM"))

    n <- nrow(fit)

    top_tags <-
        glmQLFTest(fit,coef=coef,contrast=contrast) %>%
        topTags(n=n,sort.by="none")

    pfunc <- function(i, mag) {
        if (mag == 0.0)
            top_treats <- top_tags
        else
            top_treats <-
                glmTreat(fit, coef=coef, contrast=contrast, lfc=mag) %>%
                topTags(n=n, sort.by="none")

        top_treats$table$PValue[i]
    }

    confects <- nest_confects(n, pfunc, fdr=fdr, max=max, step=step)

    logFC <- top_tags$table$logFC[confects$index]
    confects$signed_confect <- sign(logFC) * confects$confect
    confects$logFC <- logFC
    confects$logCPM <- top_tags$table$logCPM[confects$index]
    confects$name <- rownames(top_tags$table)[confects$index]

    confects
}
