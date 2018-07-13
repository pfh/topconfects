
# Ensure matrix
# Ensure column names
clean_design_matrix <- function(design) {
    design <- as.matrix(design)

    assert_that(is.numeric(design))

    if (is.null(colnames(design)))
        colnames(design) <- paste0("X", seq_len(ncol(design)))

    assert_that(length(unique(colnames(design))) == ncol(design))

    design
}


# Helper function: extract y and weights, use limma to calculate variance moderation (non-trended)
make_context <- function(object, design) {
    eawp <- limma::getEAWP(object)
    y <- eawp$exprs
    if (is.null(eawp$weights))
        weights <- matrix(1, nrow=nrow(y), ncol=ncol(y))
    else
        weights <- eawp$weights

    assert_that(identical( dim(y), dim(weights) ))

    # TODO: missing values
    assert_that(!any(is.na(y)))
    assert_that(!any(is.na(weights)))

    n_items <- nrow(y)
    n_samples <- ncol(y)
    n_coef <- ncol(design)
    assert_that(nrow(design) == n_samples)

    fit <- limma::lmFit(y, design, weights=weights) 
    fit <- limma::eBayes(fit)

    fit$genes <- eawp$probes

    list(
        y=y, 
        weights=weights, 
        design=design,
        fit=fit, 
        n_items=n_items, 
        n_samples=n_samples, 
        n_coef=n_coef)
}

context_vcov <- function(context, row) {
    s2 <- context$fit$s2.post[row]
    w <- context$weights[row,]
    X <- context$design
    Xt.W <- t(w*X)  # Multiply each row by corresponding weight, then transpose
    estimator <- solve(Xt.W %*% X) %*% Xt.W

    #Should match estimation in limma:
    #expected_coef <- as.vector(estimator %*% context$y[row,])
    #assert_that(all.equal(expected_coef, as.vector(context$fit$coefficients[row,])))

    # Variance of noise in each y
    y_var <- s2/w
    # Propagation of uncertainty
    estimator %*% (y_var*t(estimator))
}


#' Confident non-linear effect sizes using Wald tests and limma
#'
#' Fit a linear model to each gene using limma, then calculate non-linear effect sizes that are a function of the fitted coefficients.
#'
#' For RNA-Seq data, use of \code{limma::voom} is necessary to properly weight expression levels. 
#'
#' @param object An expression matrix or EList object, or anything limma's lmFit will accept.
#'
#' @param design Design matrix.
#'
#' @param effect A non-linear effect size definition, created with one of the \code{effect_...} functions.
#'
#' @param fdr False Discovery Rate to control for.
#'
#' @param step Granularity of effect sizes to test.
#'
#' @return
#'
#' See \code{\link{nest_confects}} for details of how to interpret the result.
#'
#' @export
limma_nonlinear_confects <- function(object, design, effect, fdr=0.05, step=0.001) {
    design <- clean_design_matrix(design)
    effect <- effect(colnames(design))

    # Don't currently support effect sizes that may restrict along more than one direction
    assert_that(effect$df == 1)

    context <- make_context(object, design)

    x <- numeric(context$n_items)
    se <- numeric(context$n_items)

    for(i in seq_len(context$n_items)) {
        coefs <- context$fit$coefficients[i,]
        vcovs <- context_vcov(context, i)

        result <- effect$calc(coefs)
        x[i] <- result$effect
        # Propagation of uncertainty, with local linear approximation to effect size function
        # sqrt( gradient %*% vcovs %*% t(gradient) )
        se[i] <- sqrt( sum(colSums(vcovs*result$gradient)*result$gradient) )
    }

    confects <- normal_confects(x, se, df=context$fit$df.total, fdr=fdr, step=step)

    confects$table$AveExpr <- context$fit$Amean[confects$table$index]
    if (!is.null(rownames(context$fit)))
        confects$table$name <- rownames(context$fit)[confects$table$index]
    else
        confects$table$name <- as.character(confects$table$index)

    if (!is.null(context$fit$genes)) {
        confects$table <- cbind(confects$table, context$fit$genes[confects$table$index,,drop=FALSE])
    }

    confects$limits <- effect$limits
    confects$limma_fit <- context$fit

    confects
}





