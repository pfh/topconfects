

#' Confident log2 fold changes using limma's treat function, linear effects only
#'
#' For all possible absolute log2 fold changes, which genes have at least this fold change at a specified False Discovery Rate?
#'
#' \code{fit} should be produced using \code{lmFit}, optionally \code{contrasts.fit} if a contrast is needed, and then \code{eBayes}, as in normal limma usage. 
#'
#' @param fit A limma MArrayLM object.
#'
#' @param coef Column number of coefficient or contrast to test.
#'
#' @param fdr False Discovery Rate to control for.
#'
#' @param step Granularity of log2 fold changes to test.
#'
#' @return
#'
#' See \code{\link{nest_confects}} for details of how to interpret the result.
#'
#' @export
limma_confects <- function(fit, coef=NULL, fdr=0.05, step=0.01) {
    assert_that(is(fit, "MArrayLM"), msg="fit must be an MArrayLM object")
    assert_that(!is.null(fit$df.prior), msg="Need to run eBayes first")
    
    if (is.null(coef) && ncol(fit$coefficients) == 1)
        coef <- 1

    assert_that(length(coef)==1, msg="Need to specify one coefficient to test")

    n <- nrow(fit)

    pfunc <- function(i, mag) {
        top_treats <-
            treat(fit, lfc=mag) %>%
            topTreat(coef=coef, sort.by="none",n=n)
        top_treats$P.Value[i]
    }

    confects <- nest_confects(n, pfunc, fdr=fdr, step=step)
    logFC <- fit$coefficients[confects$table$index, coef]
    confects$table$confect <- sign(logFC) * confects$table$confect
    confects$table$effect <- logFC

    confects$table$AveExpr <- fit$Amean[confects$table$index]
    if (!is.null(rownames(fit)))
        confects$table$name <- rownames(fit)[confects$table$index]
    else
        confects$table$name <- as.character(confects$table$index)
    confects$effect_desc <- "log2 fold change"

    confects$limma_fit <- fit

    if (!is.null(fit$genes)) {
        confects$table <- cbind(confects$table, fit$genes[confects$table$index,,drop=FALSE])
    }

    confects
}




