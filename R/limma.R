

#' Confident log2 fold changes based on a limma fit object
#'
#' For all possible absolute log2 fold changes (LFC), which genes have at least this fold change at a specified False Discovery Rate?
#'
#' Results are presented in a table such that for any given LFC, if the reader chooses the genes with abs(confect) less than this they are assured that this set of genes has at least this LFC (with the specified FDR). The confect column may also be viewed as a confidence bound on the LFC of each gene, with a dynamic correction for multiple testing.
#'
#' \code{fit} should be produced using \code{lmFit}, and optionally \code{contrasts.fit} if a contrast is needed. It is not necessary to use \code{eBayes}, this function calls \code{eBayes} itself.
#'
#' @param fit A limma MArrayLM object.
#'
#' @param coef Column number of coefficient or contrast to test.
#'
#' @param fdr False Discovery Rate to control for.
#'
#' @param step Granularity of log2 fold changes to test.
#' 
#' @param trend Should \code{treat(..., trend=TRUE)} be used?
#'
#' @return
#'
#' See \code{\link{nest_confects}} for details of how to interpret the result.
#'
#' @export
limma_confects <- function(fit, coef=NULL, fdr=0.05, step=0.001, trend=FALSE) {
    assert_that(is(fit, "MArrayLM"), msg="fit must be an MArrayLM object")
    
    if (is.null(coef) && ncol(fit$coefficients) == 1)
        coef <- 1

    assert_that(length(coef)==1, msg="Need to specify one coefficient to test")

    n <- nrow(fit)

    fit <- limma::eBayes(fit, trend=trend)

    # Based on limma::treat()
    #acoef <- abs(fit$coefficients[,coef])
    #se <- fit$stdev.unscaled[,coef] * sqrt(fit$s2.post)
    #df_total <- fit$df.total

    # pfunc <- function(i, mag) {
    #     acoef_i <- acoef[i]
    #     se_i <- se[i]
    #     df_total_i <- df_total[i]

    #     tstat_right <- (acoef_i-mag)/se_i
    #     tstat_left <- (acoef_i+mag)/se_i

    #     pt(tstat_right, df=df_total_i,lower.tail=FALSE) + 
    #         pt(tstat_left,df=df_total_i,lower.tail=FALSE)
    # }

    # confects <- nest_confects(n, pfunc, fdr=fdr, step=step)
    #logFC <- fit$coefficients[confects$table$index, coef]
    #confects$table$confect <- sign(logFC) * confects$table$confect
    #confects$table$effect <- logFC

    effect <- fit$coefficients[,coef]
    sd <- fit$stdev.unscaled[,coef] * sqrt(fit$s2.post)
    df <- fit$df.total
    confects <- normal_confects(effect, sd, df, fdr=fdr, step=step)

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


#
# Version of limma_confects using the limma treat function.
#
limma_confects_limma <- function(fit, coef=NULL, fdr=0.05, step=0.001, trend=FALSE) {
    assert_that(is(fit, "MArrayLM"), msg="fit must be an MArrayLM object")
    
    if (is.null(coef) && ncol(fit$coefficients) == 1)
        coef <- 1

    assert_that(length(coef)==1, msg="Need to specify one coefficient to test")

    n <- nrow(fit)

    pfunc <- function(i, mag) {
        tested <- limma::treat(fit, lfc=mag, trend=trend)
        top_treats <- limma::topTreat(tested, coef=coef, sort.by="none",n=n)
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



