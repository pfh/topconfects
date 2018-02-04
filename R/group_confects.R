
#
# Differential exon usage, 3' end points, 5' start points, etc
#
#

group_effect_1 <- function(design, coef, effect_func, design_common=NULL) {
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

# Group effect for differential splicing
#
# Create a group effect object to detect differential splicing (or similar). The effect size is the root sum of squared differences of the coefficient from the mean.
#
# The coefficient should represent a difference between two conditions.
#
# The idea is to detect differences in differential expression between features in a group (ie exons in a gene).
#
# @param design Design matrix.
#
# @param coef Column number in design matrix of the coefficient to be tested.
#
# @param design_common Experimental! Optional sample-level design matrix. For example, this can be used to account for a batch effect or matched samples.
#
# @return
#
# A group effect object.
#
# @seealso \code{\link{effect_rssm}}
#
# @export
group_effect_rssm <- function(design, coef, design_common=NULL) 
    group_effect_1(design, coef, effect_rssm, design_common=design_common)



group_effect_2 <- function(design, coef1, coef2, effect_func, design_common=NULL) {
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
        design_common = design_common,
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
# The _stepup and _stepdown versions may be used to look for shifts in start or end of transcription from RNA-seq or microarray data, where the observed levels are expected to be cumulative or reverse cumulative.
#'
#' @param design Design matrix.
#'
#' @param coef1 Column number of coefficient for first condition in design matrix.
#'
#' @param coef2 Column number of coefficient for second condition in design matrix.
#'
#' @param design_common Experimental! Optional sample-level design matrix. For example, this can be used to account for a batch effect or matched samples.
#'
#' @return
#'
#' A group effect object.
#'
#' @seealso \code{\link{effect_shift}}, \code{\link{effect_shift_log2}}
#'
#' @export
group_effect_shift <- function(design, coef1, coef2, design_common=NULL) 
    group_effect_2(design, coef1, coef2, effect_shift, design_common=design_common)

# @rdname group_effect_shift
# @export
group_effect_shift_stepup <- function(design, coef1, coef2, design_common=NULL) 
    group_effect_2(design, coef1, coef2, effect_shift_stepup, design_common=design_common)

# @rdname group_effect_shift
# @export
group_effect_shift_stepdown <- function(design, coef1, coef2, design_common=NULL) 
    group_effect_2(design, coef1, coef2, effect_shift_stepdown, design_common=design_common)

#' @rdname group_effect_shift
#' @export
group_effect_shift_log2 <- function(design, coef1, coef2, design_common=NULL) 
    group_effect_2(design, coef1, coef2, effect_shift_log2, design_common=design_common)

# @rdname group_effect_shift
# @export
group_effect_shift_stepup_log2 <- function(design, coef1, coef2, design_common=NULL) 
    group_effect_2(design, coef1, coef2, effect_shift_stepup_log2, design_common=design_common)

# @rdname group_effect_shift
# @export
group_effect_shift_stepdown_log2 <- function(design, coef1, coef2, design_common=NULL) 
    group_effect_2(design, coef1, coef2, effect_shift_stepdown_log2, design_common=design_common)


#
# Make a function to assemble a group design matrix
#
group_design_maker <- function(design, design_common=NULL) {
    stopifnot(is.matrix(design))

    if (!is.null(design_common)) {
        stopifnot(is.matrix(design_common))
        stopifnot(nrow(design_common) == nrow(design))
    }

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

        if (!is.null(design_common)) {
            design_common_rep <- do.call(rbind, rep(list(design_common), m))
            big_design <- cbind(big_design, design_common_rep)
        }

        big_design
    }

    memoise(group_design)
}

#' Group confects (differential 5' or 3' end usage, etc)
#'
#' Find differential exon usage, etc.
#'
#' If the order of the members of a group is important, ensure the order of rows is correct in the original matrix.
#'
#' Groups with less than two members will be ignored.
#'
#' To construct a group design matrix, the design from \code{fit} will be repeated in a block diagonal matrix.
#'
#' @param fit For edgeR, an edgeR DGEGLM object.
#'
#' @param object For limma, an expression matrix or EList object, or anything limma's lmFit will accept.
#'
#' @param group_id A factor of length \code{nrow(fit)}, assigning items to groups (eg genes).
#'
#' @param group_effect A group effect object created by one of the A \code{group_effect_...} functions.
#'
#' @param fdr False Discovery Rate to maintain.
#'
#' @param step Step size when calculating confident effect sizes.
#' 
#' @param trend For limma, should \code{eBayes(trend=TRUE)} be used?
#'
#' @return
#'
#' See \code{\link{nest_confects}} for details of how to interpret the result.
#'
#' @export
edger_group_confects <- function(fit, group_id, group_effect, fdr=0.05, step=0.01) {
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

    get_effect <- group_effect$get_effect
    get_effect <- memoise(get_effect)
    design <- group_effect$design
    design_common <- group_effect$design_common
    group_design <- group_design_maker(design, design_common)

    n_items <- nrow(y)
    n_samples <- ncol(y)
    n_coef <- ncol(design)
    assert_that(nrow(design) == n_samples)
    assert_that(length(offset) == n_samples)

    if (!is.null(design_common)) {
        warning("Deviance moderation is currently disabled when design_common is used.")
        fit$df.prior <- 0.0
    }

    assert_that(length(fit$df.prior) == 1)

    df_prior <- rep(fit$df.prior, n)
    # Hmm
    s2_prior_item <- broadcast(fit$var.prior,n_items)
    s2_prior <- map_dbl(members, function(indices) mean(s2_prior_item[indices]))

    df_residual <- sizes*(n_samples-n_coef)

    dispersions <- broadcast(fit$dispersion, n_items)

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





#' @rdname edger_group_confects
#' @export
limma_group_confects <- function(object, group_id, group_effect, fdr=0.05, step=0.01, trend=FALSE) {
    group_id <- factor(group_id)

    eawp <- getEAWP(object)

    y <- eawp$exprs
    if (is.null(eawp$weights))
        weights <- matrix(1, nrow=nrow(y), ncol=ncol(y))
    else
        weights <- eawp$weights

    assert_that(length(group_id) == nrow(y))

    members <- split(seq_len(nrow(y)), group_id)
    members <- members[ map_int(members,length) >= 2 ]
    sizes <- map_int(members,length)
    n <- length(members)

    get_effect <- group_effect$get_effect
    get_effect <- memoise(get_effect)
    design <- group_effect$design
    design_common <- group_effect$design_common
    group_design <- group_design_maker(design, design_common)

    if (is.null(design_common)) {
        design_full <- design
    } else {
        design_full <- cbind(design, design_common)
    }

    n_items <- nrow(y)
    n_samples <- ncol(y)
    n_coef <- ncol(design)
    assert_that(nrow(design) == n_samples)

    limma_fit <- lmFit(object, design) %>% 
        eBayes(trend=trend)

    if (!is.null(design_common)) {
        warning("Variance moderation is currently disabled when design_common is used.")
        limma_fit$df.prior <- 0.0
    }

    AveExpr <- map_dbl(members,function(items) mean(limma_fit$Amean[items]))

    assert_that(length(limma_fit$df.prior) == 1)
    df_prior <- rep(limma_fit$df.prior, n)
    # Hmm
    s2_prior_item <- broadcast(limma_fit$s2.prior, n_items)
    s2_prior <- map_dbl(members, function(indices) mean(s2_prior_item[indices]))

    df_residual <- sizes*(n_samples-n_coef)

    fit_features <- function(i, cons=NULL, equality=FALSE, initial=NULL) {
        m <- sizes[i]
        this_design <- group_design(m)
        this_y <- do.call(c, 
            lapply(members[[i]], function(j) y[j,]))
        this_weights <- do.call(c,
            lapply(members[[i]], function(j) weights[j,]))

        #(if (is.null(cons)) constrained_fit_newton else constrained_fit_slsqp)(
        constrained_fit_slsqp(
            this_y,
            this_design,
            devi_normal(this_weights),
            cons,
            initial=initial,
            equality=equality)
    }

    effect <- function(i)
        get_effect(sizes[i])

    confects <- nonlinear_confects(
        df_residual, s2_prior, df_prior, fit_features, effect, fdr, step, tech_rep=sizes)
    # Intent of tech_rep=sizes: Each *sample* only counts for one observation of the residual deviance

    confects$table$AveExpr <- AveExpr[confects$table$index]
    confects$table$name <- names(members)[confects$table$index]

    confects$limma_fit <- limma_fit
    confects$members <- members
    confects$limits <- group_effect$limits

    confects
}

# TODO: remove code duplication



