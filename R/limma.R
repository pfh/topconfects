

#' Confident log2 fold changes using limma's TREAT
#'
#' For all possible absolute log2 fold changes, which genes have at least this fold change at a specified False Discovery Rate?
#'
#' @param fit A limma MArrayLM object with one coefficient. (Use limma::contrasts.fit reduce the fit object down to the contrast or coefficient that you are interested in.)
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
limma_confects <- function(fit, fdr=0.05, max=30.0, step=0.05) {
    assert_that(is(fit, "MArrayLM"))
    assert_that(ncol(fit) == 1)

    n <- nrow(fit)

    pfunc <- function(i, mag) {
        top_treats <-
            treat(fit, lfc=mag) %>%
            topTreat(sort.by="none",n=n)
        top_treats$P.Value[i]
    }

    confects <- nest_confects(n, pfunc, fdr=fdr, max=max, step=step)
    logFC <- fit$coefficients[confects$index, 1]
    confects$signed_confect <- sign(logFC) * confects$confect
    confects$logFC <- logFC
    confects$AveExpr <- fit$Amean[confects$index]
    confects$name <- rownames(fit)[confects$index]

    confects
}




