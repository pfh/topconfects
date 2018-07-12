
#' Group effect sizes starting from limma analysis
#'
#' Fit a linear model to each feature using limma, then calculate non-linear effect sizes on groups of features that are a function of the fitted coefficients.
#'
#' For RNA-Seq data, use of \code{limma::voom} is necessary to properly weight expression levels. 
#'
#' @param object An expression matrix or EList object, or anything limma's lmFit will accept.
#'
#' @param design Design matrix.
#'
#' @param grouping A data frame with at least two columns, one named "group" and one named "name". The "name" column should correspond to row names in \code{object}. The "group" column should provide a grouping of the named rows. The order within groups may be important.
#'
#' @param effect A non-linear group effect size definition, created with one of the \code{group_effect_...} functions.
#'
#' @param intercor Correlation of the contributions of different features to the effect size. 0 implies no correlation, 1 implies perfect correlation. Larger values will produce more conservative confects. Effect size estimates are unaffected.
#'
#' @param df_conservative Use a conservative value of the df as no more than that of for single feature tests (ie df.total produced by limma). By default, a Welch-Satterthwaite approximation is used, with each member in a group contributing its down df (including prior df). This default is appropriate for gene-set enrichment, but might not be appropriate for gene effect sizes where the features are (eg) exons, and biological variation might affect the whole gene.
#'
#' @param fdr False Discovery Rate to control for.
#'
#' @param step Granularity of effect sizes to test.
#'
#' @param full Include some further statistics used to calculate confects in the output.
#'
#' @return
#'
#' See \code{\link{nest_confects}} for details of how to interpret the result.
#'
#' @export
limma_group_confects <- function(object, design, grouping, effect, 
        intercor=0.0, df_conservative=FALSE, fdr=0.05, step=0.001, full=FALSE) {
    assert_that(intercor >= 0 && intercor <= 1)

    design <- clean_design_matrix(design)
    effect <- effect(colnames(design))
    context <- make_context(object, design)

    if (is.atomic(grouping)) {
        assert_that(length(grouping) == nrow(context$n_items))
        grouping <- data.frame(
            group=grouping, 
            name=rownames(context$y), 
            stringsAsFactors=FALSE)
    }

    assert_that("group" %in% names(grouping))
    assert_that("name" %in% names(grouping))

    indices <- match(grouping$name, rownames(context$y))
    assert_that(!any(is.na(indices)))
    members <- split(indices, grouping$group)

    n_coef <- context$n_coef
    n_groups <- length(members)

    x <- numeric(n_groups)
    se <- numeric(n_groups)
    df <- numeric(n_groups)
    average <- numeric(n_groups)

    for(i in seq_len(n_groups)) {
        this_members <- members[[i]]
        this_dfs <- context$fit$df.total[this_members]
        this_coefs <- context$fit$coefficients[this_members,,drop=FALSE]
        result <- effect$calc(this_coefs)
        x[i] <- result$effect

        u2 <- numeric(length(this_members))
        for(j in seq_along(this_members)) {
            vcovs <- context_vcov(context, this_members[j])
            grads <- result$gradient[j,]
            u2[j] <- sum(colSums(vcovs*grads)*grads)
        }

        u <- sqrt(u2)
        sum_u <- sum(u)
        sum_u2 <- sum(u2)

        # If no inter-feature correlation
        # se[i] <- sqrt(sum_u2)

        # Contribution of each feature to the effect size may be correlated
        # Assume a correlation matrix with 1 on the diagonal and intercor everywhere else
        se[i] <- sqrt( sum_u2 + intercor*(sum_u^2-sum_u2) )

        # Generalized Welch-Satterthwaite after R Willink, 2007
        if (all(this_dfs == Inf))
            df[i] <- Inf
        else if (df_conservative)
            df[i] <- min(this_dfs)
        else
            df[i] <- max(
                1,
                se[i]^4 / 
                  sum( (u2+intercor*u*(sum_u-u))^2 / this_dfs )
            )

        # Average of the total expression in each sample.
        # Assumes a log2 scale!
        average[i] <- log2(sum(2^context$fit$Amean[this_members]))
    }

    # Conservative: 
    # Assumes features from the same sample subject to some sort of groupwise noise.
    # (irrespective of intercor)
    #df <- unique(context$fit$df.total)
    #assert_that(length(df) == 1, msg="This should never happen.")

    confects <- normal_confects(x, se, df, fdr=fdr, step=step, full=full)
    confects$table$AveExpr <- average[confects$table$index]
    confects$table$name <- names(members)[ confects$table$index ]

    confects$limits <- effect$limits
    confects$limma_fit <- context$fit

    confects
}



