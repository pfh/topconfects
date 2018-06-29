
#' Confident log2 fold changes based on a DESeq2 analysis
#'
#' For all possible absolute log2 fold changes, which genes have at least this fold change at a specified False Discovery Rate? This is built by repeatedly calling DESeq2::results with the "greaterAbs" alternative hypothesis. 
#'
#' Results are presented in a table such that for any given LFC, if the reader chooses the genes with abs(confect) less than this they are assured that this set of genes has at least this LFC (with the specified FDR). The confect column may also be viewed as a confidence bound on the LFC of each gene, with a dynamic correction for multiple testing.
#'
#' @param object Object produced by the \code{DESeq2::DESeq} function.
#'
#' @param ... Further arguments to \code{DESeq2::results}. At a minimum you should specify either \code{contrast=} or \code{name=}.
#'
#' @param fdr False Discovery Rate to control for.
#'
#' @param step Granularity of log2 fold changes to test.
#'
#' @return
#'
#' See \code{\link{nest_confects}} for details of how to interpret the result.
#'
#' The \code{filtered} column in the result indicates whether DESeq2 filtered the gene. Such genes do not count toward the total number of genes when controlling FDR. If your intention is to obtain a ranking of all genes, you should disable this with \code{deseq2_confects(..., cooksCutoff=Inf, independentFiltering=FALSE)}.
#'
#' @export
deseq2_confects <- function(object, ..., fdr=0.05, step=0.01) {
    n <- nrow(object)

    # Calculate effect column, determine which genes are filtered by DESeq2
    std_result <- DESeq2::results(object, ...)

    # DESeq2 may filter some genes before emitting pvalues
    filtered <- is.na(std_result$padj)
    subset_out <- which(!filtered)
    subset_n <- length(subset_out)

    pfunc <- function(subset_i, mag) {
        i <- subset_out[subset_i]

        result <- DESeq2::results(object, ..., altHypothesis="greaterAbs", lfcThreshold=mag)
        result_p <- result$pvalue[i]

        # This should never happen:
        assert_that(!any(is.na(result_p)), msg="unexpected extra filtering by DESeq2")

        result_p
    }

    confects <- nest_confects(subset_n, pfunc, fdr=fdr, step=step)
    confects$table$index <- subset_out[ confects$table$index ]

    if (subset_n < n) {
        extra <- data.frame(
            rank=(subset_n+1):n, 
            index=which(filtered), 
            confect=rep(NA,n-subset_n))
        confects$table <- rbind(confects$table, extra)
    }

    confects$table$effect <- std_result$log2FoldChange[confects$table$index]
    confects$table$confect <- sign(confects$table$effect) * confects$table$confect
    
    confects$table$baseMean <- std_result$baseMean[confects$table$index]
    confects$table$name <- rownames(object)[confects$table$index]
    confects$table$filtered <- filtered[confects$table$index]

    confects
}



