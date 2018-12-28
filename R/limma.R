

#' Confident log2 fold changes based on a limma fit object
#'
#' For all possible absolute log2 fold changes (LFC), which genes have at least
#' this fold change at a specified False Discovery Rate (FDR)?
#'
#' Results are presented in a table such that for any given LFC, if the reader
#' chooses the genes with \code{abs(confect)} less than this they are assured
#' that this set of genes has at least this LFC (with the specified FDR). Once
#' this set of genes is selected, the confect values provide confidence bounds
#' with False Coverage-statement Rate at the same level as the FDR.
#'
#' \code{fit} should be produced using \code{lmFit}. It is not necessary to use
#' \code{eBayes}, this function calls \code{eBayes} itself.
#'
#' To test contrasts, this function can also be used with the result of
#' \code{contrasts.fit}, but limma's handling of weights may be approximate (for
#' example if \code{voom} has been used). For exact results for a contrast, use
#' \code{contrastToCoef} to adjust the design matrix given to \code{lmFit}.
#'
#' @param fit A limma MArrayLM object.
#'
#' @param coef Number or name of coefficient or contrast to test.
#'
#' @param fdr False Discovery Rate to control for.
#'
#' @param step Granularity of log2 fold changes to test.
#'
#' @param trend Should \code{eBayes(trend=TRUE)} be used?
#'
#' @param full Include some further statistics used to calculate confects in the
#'   output, and also include FDR-adjusted p-values that effect size is non-zero
#'   (note that this is against the spirit of the topconfects approach).
#'
#' @return
#'
#' See \code{\link{nest_confects}} for details of how to interpret the result.
#'
#' @examples
#'
#' #Prepare a data set
#' library(NBPSeq)
#' library(edgeR)
#' library(limma)
#' data(arab)
#' dgelist <- DGEList(arab)
#' dgelist <- calcNormFactors(dgelist)
#' cpms <- cpm(dgelist, log=TRUE)
#' # Retain genes with more than a geometric mean of 2 RPM
#' # (about 5 reads per sample)
#' cpms <- cpms[rowMeans(cpms) >= 1,]
#'
#' # Fit linear model for each gene
#' treatment <- c(FALSE,FALSE,FALSE,TRUE,TRUE,TRUE)
#' batch <- factor(c(1,2,3,1,2,3))
#' design <- model.matrix(~ treatment + batch)
#' fit <- lmFit(cpms, design)
#'
#' # Calculate top confects
#' # As voom has not been used, it is necessary to use trend=TRUE
#' limma_confects(fit, "treatmentTRUE", trend=TRUE)
#'
#' @export
limma_confects <- function(
        fit, coef=NULL, fdr=0.05, step=0.001, trend=FALSE, full=FALSE) {
    assert_that(is(fit, "MArrayLM"), msg="fit must be an MArrayLM object")

    if (is.null(coef) && ncol(fit$coefficients) == 1)
        coef <- 1

    assert_that(length(coef)==1, msg="Need to specify one coefficient to test")

    n <- nrow(fit)

    fit <- limma::eBayes(fit, trend=trend)

    effect <- fit$coefficients[,coef]
    sd <- fit$stdev.unscaled[,coef] * sqrt(fit$s2.post)
    df <- fit$df.total
    confects <- normal_confects(effect, sd, df, fdr=fdr, step=step, full=full)

    confects$table$AveExpr <- fit$Amean[confects$table$index]
    if (!is.null(rownames(fit)))
        confects$table$name <- rownames(fit)[confects$table$index]
    else
        confects$table$name <- as.character(confects$table$index)
    confects$effect_desc <- "log2 fold change"

    confects$limma_fit <- fit

    if (!is.null(fit$genes)) {
        confects$table <- cbind(
            confects$table,
            fit$genes[confects$table$index,,drop=FALSE])
    }

    confects
}


#
# Version of limma_confects using the limma treat function.
#
limma_confects_limma <- function(
        fit, coef=NULL, fdr=0.05, step=0.001, trend=FALSE) {
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
        confects$table <- cbind(
            confects$table,
            fit$genes[confects$table$index,,drop=FALSE])
    }

    confects
}



