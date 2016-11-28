

#' Effect size magnitude for limma TREAT
#'
#' @param fit A limma MArrayLM object with one coefficient. (Use limma::contrasts.fit reduce the fit object down to the contrast or coefficient that you are interested in.)
#'
#' @export
limma_confects <- function(fit, fdr=0.05, max=30.0, step=0.05) {
    assert_that(is(fit, "MArrayLM"))
    assert_that(ncol(fit) == 1)

    n <- nrow(fit)

    pfunc <- function(i, mag) {
        top_treats <-
            fit %>%
            treat(lfc=mag) %>%
            topTreat(sort.by="none",n=Inf)
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




