
#
# Topconfects result class definition
#

setClass("Topconfects", representation("list"))

setMethod("show", "Topconfects", function(object) {
    cat("$table\n")
    sub <- head(object$table, 10)
    rownames(sub) <- sub$rank
    sub$rank <- NULL
    sub$index <- NULL
    print.data.frame(sub, digits=4, right=FALSE)
    if (nrow(object$table) > 10) cat("...\n")
    cat(confects_description(object))
})

first_match <- function(options, avail, default=NULL) {
    good <- options[options %in% avail]
    if (length(good) == 0) return(default)
    good[1]
}

#
# Describe some key numbers from a result (used by show()).
#
confects_description <- function(confects) {
    result <- paste0(
        sum(!is.na(confects$table$confect)),
        " of ", nrow(confects$table), " non-zero ", confects$effect_desc,
        " at FDR ", confects$fdr, "\n")

    if (!is.null(confects$df_prior)) {
        result <- paste0(result,
            "Prior df ", sprintf("%.1f", confects$df_prior), "\n")
    }

    if (!is.null(confects$edger_fit) &&
            length(confects$edger_fit$df.prior) == 1) {
        result <- paste0(result,
            "Prior df ", sprintf("%.1f", confects$edger_fit$df.prior), "\n")
    }

    if (!is.null(confects$limma_fit) &&
            length(confects$limma_fit$df.prior) == 1) {
        result <- paste0(result,
            "Prior df ", sprintf("%.1f", confects$limma_fit$df.prior), "\n")
    }

    if (!is.null(confects$edger_fit$dispersion)) {
        result <- paste0(result,
            sprintf("Dispersion %#.2g to %#.2g\n",
                min(confects$edger_fit$dispersion),
                max(confects$edger_fit$dispersion)),
            sprintf("Biological CV %.1f%% to %.1f%%\n",
                100*sqrt(min(confects$edger_fit$dispersion)),
                100*sqrt(max(confects$edger_fit$dispersion))))
    }

    result
}




#' Top confident effect sizes plot
#'
#' Create a ggplot2 object showing the confect, effect, and average expression
#' level of top features in a Topconfects object.
#'
#' For each gene, the estimated effect is shown as a dot. The confidence bound
#' is shown as a line to positive or negative infinity, showing the set of
#' non-rejected effect sizes for the feature.
#'
#' @param confects A "Topconfects" class object, as returned from
#'   limma_confects, edger_confects, etc.
#'
#' @param n Number if items to show.
#'
#' @param limits c(lower, upper) limits on x-axis.
#'
#' @return
#'
#' A ggplot2 object. Working non-interactively, you must print() this for it to
#' be displayed.
#'
#' @examples
#'
#' # Generate some random effect sizes with random accuracies
#' n <- 100
#' effect <- rnorm(n, sd=2)
#' se <- rchisq(n, df=3)^-0.5
#'
#' # Find top confident effect sizes
#' confects <- normal_confects(effect, se)
#'
#' # Plot top confident effect sizes
#' confects_plot(confects, n=30)
#'
#' @export
confects_plot <- function(confects, n=50, limits=NULL) {
    tab <- head(confects$table, n)

    mag_col <- confects$magnitude_column
    if (is.null(mag_col)) {
        mag_col <- first_match(
            c("logCPM", "AveExpr", "row_mean", "baseMean"), names(tab))
    }

    mag_desc <- confects$magnitude_desc
    if (is.null(mag_desc)) {
        mag_desc <- mag_col
    }

    name_col <- first_match(
        c("name", "index"), names(tab))

    if (identical(mag_col,"baseMean"))
        mag_scale <- "log10"
    else
        mag_scale <- "identity"

    if (is.null(limits))
        limits <- confects$limits
    
    if (is.null(limits))
        limits <- c(NA,NA)

    min_effect <- min(0, tab$effect, na.rm=TRUE)
    max_effect <- max(0, tab$effect, na.rm=TRUE)

    if (min_effect == max_effect) {
        min_effect <- -1
        max_effect <- 1
    }

    if (is.na(limits[1]) & is.na(limits[2])) {
        max_abs_effect <- max(-min_effect,max_effect)
        limits <- c(-max_abs_effect*1.05, max_abs_effect*1.05)
    } else if (is.na(limits[1])) {
        limits[1] <- min_effect * 1.05
    } else if (is.na(limits[2])) {
        limits[2] <- max_effect * 1.05
    }   

    assert_that(is.numeric(limits), length(limits) == 2)

    tab$confect_from <- limits[1]
    tab$confect_to <- limits[2]
    positive <- !is.na(tab$confect) & tab$effect > 0
    tab$confect_from[positive] <- tab$confect[positive]
    negative <- !is.na(tab$confect) & tab$effect < 0
    tab$confect_to[negative] <- tab$confect[negative]

    tab$name <- factor(tab[[name_col]],rev(tab[[name_col]]))

    p <- ggplot(tab, aes_string(y="name", x="effect")) +
        geom_vline(xintercept=0) +
        geom_segment(aes_string(
            yend="name", x="confect_from", xend="confect_to")) +
        geom_point(aes_string(size=mag_col)) +
        scale_x_continuous(expand=c(0,0), limits=limits, oob=function(a,b) a) +
        labs(x = confects$effect_desc, y="", size=mag_desc) +
        theme_bw()

    if (identical(mag_col,"baseMean"))
        p <- p + scale_size(trans="log10")

    p
}

#' Mean-expression vs effect size plot (deprecated)
#'
#' Note: I now recommend using \code{plot_confects_me2} instead of this plot. 
#' Like plotMD in limma, plots effect size against mean expression level.
#' However shows "confect" on the y axis rather than "effect" ("effect" is shown
#' underneath in grey). This may be useful for assessing whether effects are
#' only being detected in highly expressed genes.
#'
#' @param confects A "Topconfects" class object, as returned from
#'   \code{limma_confects}, \code{edger_confects}, or \code{deseq2_confects}.
#'
#' @return
#'
#' The two types of points in this plot make it quite confusing to explain.
#' \code{plot_confects_me2} is recommended instead.
#'
#' A ggplot2 object. Working non-interactively, you must print() this for it to
#' be displayed.
#'
#' @examples
#'
#' library(NBPSeq)
#' library(edgeR)
#' library(limma)
#'
#' data(arab)
#'
#' # Extract experimental design from sample names
#' treat <- factor(substring(colnames(arab),1,4), levels=c("mock","hrcc"))
#' time <- factor(substring(colnames(arab),5,5))
#'
#' # Keep genes with at least 3 samples having an RPM of more than 2
#' y <- DGEList(arab)
#' keep <- rowSums(cpm(y)>2) >= 3
#' y <- y[keep,,keep.lib.sizes=FALSE]
#' y <- calcNormFactors(y)
#'
#' # Find top confident fold changes by topconfects-limma-voom method
#' design <- model.matrix(~time+treat)
#' voomed <- voom(y, design)
#' fit <- lmFit(voomed, design)
#' confects <- limma_confects(fit, "treathrcc")
#'
#' # Plot confident effect size against mean expression
#' # (estimated effect size also shown as grey dots)
#' confects_plot_me(confects)
#'
#' @export
confects_plot_me <- function(confects) {
    tab <- confects$table

    mag_col <- confects$magnitude_column
    if (is.null(mag_col)) {
        mag_col <- first_match(
            c("logCPM", "AveExpr", "row_mean", "baseMean"), names(tab))
    }

    assert_that(!is.null(mag_col), msg="No mean expression column available.")

    mag_desc <- confects$magnitude_desc
    if (is.null(mag_desc)) {
        mag_desc <- mag_col
    }
    
    non_na_tab <- tab[!is.na(tab$confect),]
    non_na_tab$color <- ifelse(non_na_tab$effect >= 0, "#cc0000", "#0000cc")
    
    p <- ggplot(tab, aes_string(x=mag_col)) +
        geom_point(aes_string(y="effect"), color="#cccccc") +
        geom_hline(yintercept=0) +
        geom_point(data=non_na_tab, mapping=aes_string(y="confect"), color=non_na_tab$color) +
        labs(x=mag_desc, y=confects$effect_desc) +
        theme_bw()

    if (identical(mag_col,"baseMean"))
        p <- p + scale_x_continuous(trans="log10")

    p
}



#' Mean-expression vs effect size plot (version 2)
#'
#' Like plotMD in limma, plots effect size against mean expression level.
#' Confect values are indicated using color.
#' This may be useful for assessing whether effects are
#' only being detected in highly expressed genes.
#'
#' @param confects A "Topconfects" class object, as returned from
#'   \code{limma_confects}, \code{edger_confects}, or \code{deseq2_confects}.
#'
#' @param breaks A vector of confect thresholds to color points with.
#'   Chooses a sensible set of breaks if none are given.
#'
#' @return
#'
#' A ggplot2 object. Working non-interactively, you must print() this for it to
#' be displayed.
#'
#' @examples
#'
#' library(NBPSeq)
#' library(edgeR)
#' library(limma)
#'
#' data(arab)
#'
#' # Extract experimental design from sample names
#' treat <- factor(substring(colnames(arab),1,4), levels=c("mock","hrcc"))
#' time <- factor(substring(colnames(arab),5,5))
#'
#' # Keep genes with at least 3 samples having an RPM of more than 2
#' y <- DGEList(arab)
#' keep <- rowSums(cpm(y)>2) >= 3
#' y <- y[keep,,keep.lib.sizes=FALSE]
#' y <- calcNormFactors(y)
#'
#' # Find top confident fold changes by topconfects-limma-voom method
#' design <- model.matrix(~time+treat)
#' voomed <- voom(y, design)
#' fit <- lmFit(voomed, design)
#' confects <- limma_confects(fit, "treathrcc")
#'
#' # Plot confident effect size against mean expression
#' # (estimated effect size also shown as grey dots)
#' confects_plot_me2(confects)
#'
#' @export
confects_plot_me2 <- function(confects, breaks=NULL) {
    tab <- confects$table
    
    mag_col <- confects$magnitude_column
    if (is.null(mag_col)) {
        mag_col <- topconfects:::first_match(
            c("logCPM", "AveExpr", "row_mean", "baseMean"), names(tab))
    }
    
    assert_that(!is.null(mag_col), msg="No mean expression column available.")
    
    mag_desc <- confects$magnitude_desc
    if (is.null(mag_desc)) {
        mag_desc <- mag_col
    }
    
    # Provide sensible breaks if not given
    if (length(breaks) == 0) {
        breaks <- 0.0
        max_confect <- max(0.0, abs(tab$confect), na.rm=TRUE)
        big_step <- 0.1
        repeat {
            step <- big_step
            if (step*4 >= max_confect) break
            step <- big_step*2
            if (step*4 >= max_confect) break
            step <- big_step*5
            if (step*4 >= max_confect) break
            big_step <- big_step*10
        }
        breaks <- seq(0, max_confect, by=step)
    }
    
    tab$color <- rep(NA, nrow(tab))
    for(item in sort(breaks)) {
        tab$color[ !is.na(tab$confect) & abs(tab$confect) >= item ] <- item
    }
    tab$color <- factor(tab$color, rev(breaks))
    levels(tab$color) <- paste0("> ", levels(tab$color))
    
    # Reverse order
    tab <- purrr::map_df(tab, rev)
    
    p <- ggplot(tab, aes_string(x=mag_col, y="effect", color="color")) +
        geom_hline(yintercept=0) +
        geom_point(show.legend=TRUE) +
        scale_color_discrete(
            breaks=levels(tab$color),
            na.value="#bbbbbb",
            h=c(180,380), direction=-1, drop=FALSE) +
        guides(color=guide_legend(override.aes=list(size=4))) +
        labs(
            x=mag_desc, 
            y=confects$effect_desc,
            color=paste0("|",confects$effect_desc,"|\n\nconfidently")) +
        theme_bw()
    
    if (identical(mag_col,"baseMean"))
        p <- p + scale_x_continuous(trans="log10")
    
    p
}



#' A plot to compare two rankings
#'
#' This is useful, for example, when comparing different methods of ranking
#' potentially interesting differentially expressed genes.
#'
#' @param vec1 A vector of names.
#'
#' @param vec2 Another vector of names.
#'
#' @param label1 A label to go along with vec1.
#'
#' @param label2 A label to go along with vec2.
#'
#' @param n Show at most the first n names in vec1 and vec2.
#'
#' @return
#'
#' A ggplot2 object. Working non-interactively, you must print() this for it to
#' be displayed.
#'
#' @examples
#'
#' a <- sample(letters)
#' b <- sample(letters)
#' rank_rank_plot(a,b, n=20)
#'
#' @export
rank_rank_plot <- function(
        vec1, vec2, label1="First ranking", label2="Second ranking", n=40) {
    vec1 <- as.character( head(vec1, n) )
    vec2 <- as.character( head(vec2, n) )

    df1 <- data.frame(
        rank1=seq_len(length(vec1)), name=vec1, stringsAsFactors=FALSE)
    df2 <- data.frame(
        rank2=seq_len(length(vec2)), name=vec2, stringsAsFactors=FALSE)
    link <- merge(df1, df2, by="name")
    p <- ggplot(link) +
        geom_segment(aes_string(x="1", xend="2", y="rank1", yend="rank2")) +
        geom_point(aes_string(x="1", y="rank1")) +
        geom_point(aes_string(x="2", y="rank2")) +
        geom_text(data=df1, aes_string(x="0.9",y="rank1",label="name"),
            hjust=1,vjust=0.5) +
        geom_text(data=df2, aes_string(x="2.1",y="rank2",label="name"),
            hjust=0,vjust=0.5) +
        scale_x_continuous(limits=c(0,3), breaks=c(1,2),
            minor_breaks=NULL, labels=c(label1,label2)) +
        scale_y_continuous(breaks=seq_len(n), minor_breaks=NULL,
            trans="reverse") +
        labs(x="",y="") +
        theme_bw()

    p
}

