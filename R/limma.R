

#' Confident log2 fold changes using limma's treat function, linear effects only
#'
#' For all possible absolute log2 fold changes, which genes have at least this fold change at a specified False Discovery Rate?
#'
#' @param fit A limma MArrayLM object.
#'
#' @param coef Coefficient to test, as per \code{glmTreat}. Use either coef or contrast, or have fit be the result of contrasts.fit.
#'
#' @param contrast Contrast to test, as per \code{glmTreat}. Use either coef or contrast, or have fit be the result of contrasts.fit.
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
limma_confects <- function(fit, coef=NULL, contrast=NULL, fdr=0.05, max=30.0, step=0.01) {
    assert_that(is(fit, "MArrayLM"))

    cfit <- fit
    if (!is.null(coef) || !is.null(contrast))
        cfit <- contrasts.fit(cfit, contrasts=contrast, coefficients=coef) %>%
            eBayes()

    if (is.null(fit$df.prior))
        cfit <- eBayes(cfit)

    assert_that(ncol(cfit) == 1)

    n <- nrow(cfit)

    pfunc <- function(i, mag) {
        top_treats <-
            treat(cfit, lfc=mag) %>%
            topTreat(sort.by="none",n=n)
        top_treats$P.Value[i]
    }

    confects <- nest_confects(n, pfunc, fdr=fdr, max=max, step=step)
    logFC <- cfit$coefficients[confects$table$index, 1]
    confects$table$confect <- sign(logFC) * confects$table$confect
    confects$table$effect <- logFC

    confects$table$AveExpr <- cfit$Amean[confects$table$index]
    confects$table$name <- rownames(cfit)[confects$table$index]
    confects$effect_desc <- "log2 fold change"

    confects$limma_fit <- fit

    if (!is.null(fit$genes)) {
        confects$table <- cbind(confects$table, fit$genes[confects$table$index,])
    }

    confects
}




