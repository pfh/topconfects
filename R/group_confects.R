
#
# Differential exon usage, 3' end points, 5' start points, etc
#
#

group_effect_1 <- function(design, coef, effect_func) {
    assert_that(length(coef) == 1, coef >= 1, coef <= ncol(design)) 

    get_effect <- function(n) {
        s0 <- (seq_len(n)-1)*ncol(design)
        group_coef <- coef + s0
        effect_func(group_coef)
    }

    list(
        design = design,
        signed = get_effect(2)$signed,
        limits = get_effect(2)$limits,
        get_effect = get_effect)
}

#' Group effect for differential splicing
#'
#' Create a group effect object to detect differential splicing (or similar). The effect size is the root sum of squared differences of the coefficient from the mean.
#'
#' The coefficient should represent a difference between two conditions.
#'
#' The idea is to detect differences in differential expression between features in a group (ie exons in a gene).
#'
#' @param design Design matrix.
#'
#' @param coef Column number in design matrix of the coefficient to be tested.
#'
#' @return
#'
#' A group effect object.
#'
#' @seealso \code{\link{effect_rssm}}
#'
#' @export
group_effect_rssm <- function(design, coef) group_effect_1(design, coef, effect_rssm)



group_effect_2 <- function(design, coef1, coef2, effect_func) {
    assert_that(length(coef1) == 1, coef1 >= 1, coef1 <= ncol(design)) 
    assert_that(length(coef2) == 1, coef2 >= 1, coef2 <= ncol(design))

    get_effect <- function(n) {
        s0 <- (seq_len(n)-1)*ncol(design)
        group_coef1 <- coef1 + s0
        group_coef2 <- coef2 + s0
        effect_func(group_coef1, group_coef2)
    }

    list(
        design = design,
        signed = get_effect(2)$signed,
        limits = get_effect(2)$limits,
        get_effect = get_effect)
}

#' Group effect for shifts in 5' or 3' end usage
#'
#' Create a group effect object to detect a shift in usage of features which are ordered within a group.
#'
#' group_effect_shift_log2 is almost certainly the version you want.
#'
#' The coefficients should represent the expression levels in two different conditions.
#'
#' @param design Design matrix.
#'
#' @param coef1 Column number of coefficient for first condition in design matrix.
#'
#' @param coef2 Column number of coefficient for second condition in design matrix.
#'
#' @return
#'
#' A group effect object.
#'
#' @seealso \code{\link{effect_shift}}, \code{\link{effect_shift_log2}}
#'
#' @export
group_effect_shift <- function(design, coef1, coef2) group_effect_2(design, coef1, coef2, effect_shift)

#' @rdname group_effect_shift
#' @export
group_effect_shift_log2 <- function(design, coef1, coef2) group_effect_2(design, coef1, coef2, effect_shift_log2)



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
#' @param group_effect A group effect object created by one of the A \code{group_effect_...} functions.
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

    if (is.function(group_effect))
        # Deprecated
        get_effect <- group_effect
    else {
        get_effect <- group_effect$get_effect
        design <- group_effect$design
    }
    get_effect <- memoise(get_effect)

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
        get_effect(sizes[i])

    confects <- nonlinear_confects(
        df_residual, s2_prior, df_prior, fit_features, effect, fdr, step, tech_rep=sizes)
    # Intent of tech_rep=sizes: Each *sample* only counts for one observation of the residual deviance

    confects$table$logCPM <- logCPM[confects$table$index]
    confects$table$name <- names(members)[confects$table$index]

    confects$edger_fit <- fit
    confects$members <- members
    if (is.list(group_effect))
        confects$limits <- group_effect$limits

    confects
}








