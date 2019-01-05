
# TODO: there are valid upper bounds for this as well

# TODO: proportion weighter


p_in_ball_function <- function(x,S,df, tol=1e-7) {
    # See http://logarithmic.net/pfh-files/blog/01539390477/ballhyp.html

    # Remove missing values (ie treat them as zero)
    present <- !is.na(x)
    x <- x[present]
    S <- S[present,present,drop=FALSE]
    if (length(x) == 0)
        return(function(radius) 1)

    # Rotate using eigendecomposition to point of view with no covariance
    decomp <- eigen(S, symmetric=TRUE)
    x <- t(decomp$vectors) %*% x
    d <- decomp$values
    # d is the variance of each element of x

    # Allow for redundant contrasts by not counting empty dimensions as degrees of freedom
    d_floor <- tol*max(d)
    df1 <- sum(d >= d_floor)
    d <- pmax(d, d_floor)

    x2 <- sum(x*x)
    
    pfunc <- function(radius) {
        r2 <- radius*radius
        if (x2 <= r2) return(1.0)

        # Choose w to produce x0 with desired radius
        score <- function(w) sum((x*(w/(1-w+w*d)))^2)-r2
        w <- uniroot(score, c(0,1))$root

        # Worst case null-hypothesis point
        x0 <- x*(w/(1-w+w*d))

        chisq <- sum((x-x0)^2/d)
        pf(chisq/df1, df1, df, lower.tail=FALSE)
    }

    list(df1=df1, pfunc=pfunc)
}

#' 
#'
#' @export
ball_confects <- function(effect, covar, df=Inf, fdr=0.05, step=0.001, full=FALSE) {
    n <- length(effect)

    # Sanity check
    assert_that(is.list(effect))
    assert_that(is.list(covar))
    assert_that(length(covar) == n)
    for(i in seq_len(n)) {
        assert_that(nrow(covar[[i]]) == length(effect[[i]]))
        assert_that(ncol(covar[[i]]) == length(effect[[i]]))
    }
    df <- broadcast(df, n)

    pfuncs <- map(seq_len(n), function(i) {
        p_in_ball_function(effect[[i]], covar[[i]], df[i])
    })

    pfunc <- function(indices, mag) {
        map_dbl(indices, function(i) {
            pfuncs[[i]]$pfunc(mag)
        })
    }

    confects <- nest_confects(n, pfunc, fdr=fdr, step=step, full=full)
    
    confects$table$effect <- map_dbl(
        effect[confects$table$index], ~sqrt(sum(.*.)))
    confects$table$estimates <- effect[confects$table$index]

    # FCR correction assuming number of discoveries is number 
    # of items with confect >= the current item's confect
    # Note: confidence intervals overlapping zero are counted as confidence intervals
    #       so this may be somewhat liberal if these are then ignored
    # Note: we don't do anything about different items having different sizes, so this is FCR
    #       control assuming we uniformly randomly pick a discovery 
    #       and then uniformly randomly pick an interval
    n_discoveries <- rank(-confects$table$confect, na.last=TRUE, ties.method="max")
    lower <- rep(list(NULL), n)
    upper <- rep(list(NULL), n)
    for(i in seq_len(n)) {
        j <- confects$table$index[i]
        if (is.na(confects$table$confect[i])) {
            width <- rep(NA_real_, length(effect[[j]]))
        } else {
            alpha <- fdr * n_discoveries[i]/n
            width <- sqrt(diag(covar[[j]])) * 
                qt(0.5*alpha, df[j], lower.tail=FALSE)
        }
        lower[[i]] <- effect[[j]] - width
        upper[[i]] <- effect[[j]] + width
    }

    confects$table$lower <- lower
    confects$table$upper <- upper

    if (full) {
        confects$table$covar <- covar[confects$table$index]
        confects$table$df1 <- map_dbl(pfuncs, "df1")[confects$table$index]
        confects$table$df2 <- df[confects$table$index]
        fdr_zero <- confects$table$fdr_zero
        confects$table$fdr_zero <- NULL
        confects$table$fdr_zero <- fdr_zero
    }
    
    confects$limits <- c(0,NA)

    confects
}



weighted_fit_and_contrast <- function(X,y,w,contrasts) {
    present <- (w > 0) & !is.na(y) & !is.na(w)

    sw <- sqrt(w[present])
    X <- X[present,,drop=F] * sw
    y <- y[present] * sw
    # Problem is now reduced to OLS

    #estimator <- tcrossprod(solve(crossprod(X,X)),X)
    estimator <- qr.solve(X, diag(nrow(X)))
    # will produce error on singular matrix (with default tol=1e-7)

    coef <- as.vector(estimator %*% y)
    residuals <- y - as.vector(X %*% coef)

    contrast_estimator <- crossprod(contrasts, estimator)
    effect <- as.vector(contrast_estimator %*% y)
    covar_unscaled <- tcrossprod(contrast_estimator, contrast_estimator)

    list(
        effect=effect,
        covar_unscaled=covar_unscaled,
        df=length(y) - length(coef),
        rss=sum(residuals*residuals)
    )
}

#' @export
lm_ball_confects <- function(
        y, design, contrasts, weights=NULL, 
        squeeze_var=TRUE, fdr=0.05, step=0.001, full=FALSE) {
    eawp <- limma::getEAWP(y)
    y <- eawp$exprs
    infos <- eawp$probes
    if (is.null(weights)) {
        if (is.null(eawp$weights)) {
            weights <- matrix(1, nrow=nrow(y), ncol=ncol(y))
        } else {
            weights <- eawp$weights
        }
    } else {
        if (!is.null(eawp$weights)) {
            warning("Weights in \"y\" overridden by parameter \"weights\"")
        }
    }
    
    assert_that(is.matrix(y))
    assert_that(is.matrix(weights))
    assert_that(is.matrix(design))
    assert_that(is.matrix(contrasts))
    n <- nrow(y)
    m <- ncol(y)
    assert_that(nrow(weights) == n)
    assert_that(ncol(weights) == m)
    assert_that(nrow(design) == m)
    p <- ncol(design)
    assert_that(nrow(contrasts) == p)

    fits <- map(seq_len(n), function(i) {
        weighted_fit_and_contrast(design, y[i,], weights[i,], contrasts)
    })

    df <- map_dbl(fits, "df")
    rss <- map_dbl(fits, "rss")

    scale2 <- rss/df

    if (squeeze_var) {
        squeeze <- limma::squeezeVar(scale2, df)
        df <- df + squeeze$df.prior
        scale2 <- squeeze$var.post 
    }

    effect <- map(fits, "effect")
    covar <- map(seq_len(n), function(i) {
        fits[[i]]$covar_unscaled * scale2[i]
    })

    confects <- ball_confects(
        effect, covar, df, fdr=fdr, step=step, full=full)

    confects$table$average <- map_dbl(confects$table$index, function(i) {
        mean(y[i,weights[i,]>0], na.rm=TRUE)
    })
    if (!is.null(rownames(y))) {
        confects$table$name <- rownames(y)[confects$table$index]
    } else {
        confects$table$name <- as.character(confects$table$index)
    }

    if (!is.null(infos)) {
        confects$table <- cbind(
            confects$table, infos[confects$table$index,,drop=FALSE])
    }

    if (squeeze_var) {
        confects$squeeze <- squeeze
    }
    confects$design <- design
    confects$contrasts <- contrasts

    confects
}





