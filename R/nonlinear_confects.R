
# Not vectorized!
moderated_pf <- function(f, df, s2_residual, df_residual, s2_prior, df_prior) {
    if (df_prior == Inf)
         return( pchisq(df*f/s2_prior, df=df, lower.tail=FALSE) )

    df_post <- df_prior + df_residual
    s2_post <- (s2_prior*df_prior+s2_residual*df_residual)/df_post

    if (s2_post == 0) return( 1.0 )

    return( pf(f/s2_post, df1=df, df2=df_post, lower.tail=FALSE) )
}


broadcast <- function(vec, n) {
    names(vec) <- NULL
    if (length(vec) == 1) return(rep(vec,n))
    assert_that(length(vec) == n)
    vec
}


# Called from edger_confects
edger_nonlinear_confects <- function(data, effect, fdr=0.05, max=30.0, step=0.05) {
    assert_that(is(data, "DGEGLM"))

    # Mimic edgeR's shrunk coefficients
    y <- addPriorCount(data$counts, offset=data$offset, prior.count=0.125)$y
    offset <- c(data$offset) / log(2)
    design <- data$design

    n_items <- nrow(y)
    n_samples <- ncol(y)
    n_coef <- ncol(design)
    assert_that(nrow(design) == n_samples)
    assert_that(length(offset) == n_samples)

    dispersions <- broadcast(data$dispersion, n_items)

    if (!is.null(data$df.prior)) {
        #glmQLFit
        s2_prior <- broadcast(data$var.prior, n_items)
        df_prior <- broadcast(data$df.prior, n_items)
    } else {
        #glmFit
        s2_prior <- broadcast(1.0, n_items)
        df_prior <- broadcast(Inf, n_items)
    }

    df_residual <- rep(n_samples - n_coef, n_items)

    fit <- function(i, cons=NULL, initial=NULL) {
        (if (is.null(cons)) constrained_fit_newton else constrained_fit_slsqp)(
            y[i,],
            design,
            devi_link_log2(devi_nbinom(rep(dispersions[i], n_samples))),
            cons,
            offset=offset, initial=initial)
    }

    confects <- nonlinear_confects(df_residual, s2_prior, df_prior, fit, effect, fdr, max, step)

    confects$table$logCPM <- data$AveLogCPM[confects$table$index]
    confects$table$name <- rownames(data)[confects$table$index]

    confects
}


#' Confident non-linear effect sizes using limma
#'
#' @param object An expression matrix or EList object, or anything limma's lmFit will accept.
#'
#' @param design Design matrix.
#'
#' @param effect A non-linear effect, created with one of the effect_... functions.
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
limma_nonlinear_confects <- function(object, design, effect, fdr=0.05, max=30.0, step=0.05) {
    # TODO: missing values
    eawp <- getEAWP(object)
    limma_fit <- lmFit(object, design) %>% eBayes()

    y <- eawp$exprs
    if (is.null(eawp$weights))
        weights <- matrix(1, nrow=nrow(y), ncol=ncol(y))
    else
        weights <- eawp$weights

    n_items <- nrow(y)
    n_samples <- ncol(y)
    n_coef <- ncol(design)
    assert_that(nrow(design) == n_samples)

    df_prior <- broadcast(limma_fit$df.prior, n_items)
    s2_prior <- broadcast(limma_fit$s2.prior, n_items)

    df_residual <- rep(n_samples - n_coef, n_items)

    fit <- function(i, cons=NULL, initial=NULL) {
        (if (is.null(cons)) constrained_fit_newton else constrained_fit_slsqp)(
            y[i,],
            design,
            devi_normal(weights[i,]),
            cons,
            initial=initial)
    }

    confects <- nonlinear_confects(df_residual, s2_prior, df_prior, fit, effect, fdr, max, step)

    confects$table$AveExpr <- limma_fit$Amean[confects$table$index]
    confects$table$name <- rownames(limma_fit)[confects$table$index]

    confects$limma_fit <- limma_fit

    confects
}


# Workhorse function
nonlinear_confects <- function(df_residual, s2_prior, df_prior, fit, effect, fdr, max, step) {
    n_items <- length(df_residual)

    h1_fits <- lapply(seq_len(n_items), fit)

    effects <- map_dbl(h1_fits, function(item) effect$calc( item$beta ))

    #hneg1_constraint <- effect$constraint(0)
    #hneg1_fits <- lapply(seq_len(n_items), function(i) fit(i, hneg1_constraint, h1_fits[[i]]$beta))

    pfunc <- function(indices, effect_size) {
        pos_constraint <- effect$constraint(effect_size)
        if (effect_size != 0 && effect$signed)
            neg_constraint <- effect$constraint(-effect_size)

        p <- rep(1, length(indices))
        for(i in seq_along(indices)) {
            j <- indices[i]

            h1_fit <- h1_fits[[j]]

            if (effect_size == 0 || !effect$signed) {
                # p=1 if ML estimate lies within H0 set
                if (abs(effects[j]) < effect_size) next

                h0_fit <- fit(j, pos_constraint, h1_fit$beta)

                p[i] <- moderated_pf(
                    f = (h0_fit$deviance-h1_fit$deviance)/effect$df,
                    df = effect$df,
                    s2_residual = h1_fit$deviance / df_residual[j],
                    df_residual = df_residual[j],
                    s2_prior = s2_prior[j],
                    df_prior = df_prior[j])
            } else {
                # TREAT-style p values
                # Can be < 1 even when ML estimate lies within H0 set
                h0_fit_neg <- fit(j, neg_constraint, h1_fit$beta)
                h0_fit_pos <- fit(j, pos_constraint, h1_fit$beta)
                p[i] <- 0.5*(
                    moderated_pf(
                        f = (h0_fit_neg$deviance-h1_fit$deviance)/effect$df,
                        df = effect$df,
                        s2_residual = h1_fit$deviance / df_residual[j],
                        df_residual = df_residual[j],
                        s2_prior = s2_prior[j],
                        df_prior = df_prior[j]) +
                    moderated_pf(
                        f = (h0_fit_pos$deviance-h1_fit$deviance)/effect$df,
                        df = effect$df,
                        s2_residual = h1_fit$deviance / df_residual[j],
                        df_residual = df_residual[j],
                        s2_prior = s2_prior[j],
                        df_prior = df_prior[j]))

                if (abs(effects[j]) < effect_size)
                    p[i] <- 1 - p[i]
            }
        }

        p
    }

    confects <- nest_confects(n_items, pfunc, fdr=fdr, max=max, step=step)

    if (effect$signed)
        confects$table$confect <- sign(effects[confects$table$index]) * confects$table$confect
    confects$table$effect <- effects[confects$table$index]
    confects$h1_fits <- h1_fits

    confects
}

