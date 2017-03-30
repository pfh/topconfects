
#
# Differential exon usage, 3' end points, 5' start points, etc
#
#



#' Group confects (differential exon usage, etc)
#'
#' Find differential exon usage, etc.
#'
#' If the order of the members of a group is important, ensure the order of rows is correct in the original matrix.
#'
#' Groups with less than two members will be ignored.
#'
#' To construct a group design matrix, the design from \code{fit} will be repeated in a block diagonal matrix.
#'
#' @param fit An edgeR DGEGLM object.
#'
#' @param group_id A factor of length \code{nrow(fit)}, assigning items to groups (eg genes).
#'
#' @param group_effect A \code{function(n)} returning an effect size object for a group with \code{n} members.
#'
#' @param fdr False Discovery Rate to maintain.
#'
#' @param step Step size when calculating confident effect sizes.
#'
#' @param design The design matrix is usually taken from \code{fit}, however it may be overridden with this parameter. Note that dispersion estimates and df.prior are taken from \code{fit}. So you may estimate dispersions with a simpler model than is specified here, in order to at least obtain pessimistic confects when you lack replicates.
#' 
#' @return
#'
#' See \code{\link{nest_confects}} for details of how to interpret the result.
#'
#' @export
edger_group_confects <- function(fit, group_id, group_effect, fdr=0.05, step=0.01, design=NULL) {
    group_id <- factor(group_id)
    group_effect <- memoise(group_effect)

    assert_that(is(fit, "DGEGLM"))
    assert_that(length(group_id) == nrow(fit))

    members <- split(seq_len(nrow(fit)), group_id)
    members <- members[ map_int(members,length) >= 2 ]
    sizes <- map_int(members,length)
    n <- length(members)

    logCPM <- map_dbl(members,function(items) log2(sum(2^fit$AveLogCPM[items])))

    # Mimic edgeR's shrunk coefficients
    y <- addPriorCount(fit$counts, offset=fit$offset, prior.count=0.125)$y
    offset <- c(fit$offset) / log(2)
    if (is.null(design))
        design <- fit$design

    n_items <- nrow(y)
    n_samples <- ncol(y)
    n_coef <- ncol(design)
    assert_that(nrow(design) == n_samples)
    assert_that(length(offset) == n_samples)

    assert_that(length(fit$df.prior) == 1)

    df_prior <- rep(fit$df.prior, n)
    # Hmm
    s2_prior_item <- broadcast(fit$var.prior,n_items)
    s2_prior <- map_dbl(members, function(indices) mean(s2_prior_item[indices]))

    df_residual <- sizes*(n_samples-n_coef)

    dispersions <- broadcast(fit$dispersion, n_items)

    group_design <- function(m) {
        zero_block <- matrix(0,nrow=nrow(design),ncol=ncol(design))
        big_design <- matrix(0,nrow=0,ncol=m*ncol(design))
        for(i in seq_len(m)) {
            big_design <- rbind(big_design,
                do.call(cbind,c(
                    rep(list(zero_block), i-1),
                    list(design),
                    rep(list(zero_block), m-i)
                ))
            )
        }

        big_design
    }
    group_design <- memoise(group_design)

    fit_features <- function(i, cons=NULL, equality=FALSE, initial=NULL) {
        m <- sizes[i]
        this_design <- group_design(m)
        this_y <- do.call(c, 
            lapply(members[[i]], function(j) y[j,]))
        this_offset <- rep(offset, m)
        this_dispersions <- rep(
            dispersions[members[[i]]],
            each=ncol(y))

        #(if (is.null(cons)) constrained_fit_newton else constrained_fit_slsqp)(
        constrained_fit_slsqp(
            this_y,
            this_design,
            devi_link_log2(devi_nbinom(this_dispersions)),
            cons,
            offset=this_offset, initial=initial,
            equality=equality)
    }

    effect <- function(i)
        group_effect(sizes[i])

    confects <- nonlinear_confects(
        df_residual, s2_prior, df_prior, fit_features, effect, fdr, step, tech_rep=sizes)
    # Intent of tech_rep=sizes: Each *sample* only counts for one observation of the residual deviance

    confects$table$logCPM <- logCPM[confects$table$index]
    confects$table$name <- names(members)[confects$table$index]

    confects$edger_fit <- fit
    confects$members <- members

    confects
}








