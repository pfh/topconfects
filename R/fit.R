

devi_normal <- function(weights) {
    weights2 <- 2*weights
    list(
        y_to_x = identity,
        devi = function(x,y) {
            xmy <- x-y
            weights*(xmy*xmy)
        },
        devi1 = function(x,y) {
            xmy <- x-y
            list(
                devi =  weights*(xmy*xmy),
                devi_x = weights2*xmy
            )
        },
        devi2 = function(x,y) {
            xmy <- x-y
            list(
                devi =  weights*(xmy*xmy),
                devi_x = weights2*xmy,
                devi_xx = weights2
            )
        }
    )
}


# Deviance where the variance is x+d*x^2, everything below derives from this.
# This includes the Negative Binomial distribution (but same results would be obtained if y is considered continuous)
#
devi_nbinom <- function(dispersions) {
    d <- dispersions

    list(
        y_to_x = identity,
        devi = function(x,y) {
            (-2*y*log(x)+2*(y*d+1)*log(d*x+1)/d
             +2*y*log(y)-2*(y*d+1)*log(d*y+1)/d)
        },
        devi1 = function(x,y) {
            list(
                devi = (-2*y*log(x)+2*(y*d+1)*log(d*x+1)/d
                        +2*y*log(y)-2*(y*d+1)*log(d*y+1)/d),
                devi_x = -2 * (y-x) / (x+d*x*x)
            )
        },
        devi2 = function(x,y) {
            list(
                devi = (-2*y*log(x)+2*(y*d+1)*log(d*x+1)/d
                        +2*y*log(y)-2*(y*d+1)*log(d*y+1)/d),
                devi_x = -2 * (y-x) / (x+d*x*x),
                devi_xx = -2 * (d*x*x-2*d*x*y-y) / (d*x*x+x)**2
            )
        }
    )
}


#TODO: generic quadratic variance
#TODO: tail-length


ln2 <- log(2)

devi_link_log2 <- function(devi) {
    list(
        devi = function(x,y)
            devi$devi(exp(x*ln2), y),

        devi1 = function(x,y) {
            e <- exp(x*ln2)
            out <- devi$devi1(e,y)
            list(
                devi = out$devi,
                devi_x = ln2*e*out$devi_x
            )
        },

        devi2 = function(x,y) {
            e <- exp(x*ln2)
            out <- devi$devi2(e,y)
            list(
                devi = out$devi,
                devi_x = ln2*e*out$devi_x,
                devi_xx = ln2*ln2*e*( out$devi_x + e*out$devi_xx )
            )
        },

        y_to_x = function(y) log2(devi$y_to_x(y))
    )
}


#' Simple linear contrast effect.
#'
#' This is included for completeness and testing. limma_confects() and edger_confects() provide a faster version of this using limma's treat and edgeR's topTreat functions.
#'
#' @param contrast A vector of the same length as the number of coefficients. The effect size is the dot product of this with the coefficients.
#'
#' @return
#'
#' An object defining how to calculate an effect size.
#'
#' @export
effect_contrast <- function(contrast) {
    list(
        signed = TRUE,
        df = 1,

        calc = function(beta) sum(beta * contrast),

        constraint = function(effect_size) {
            function(beta) list(
                score = sum(beta*contrast)-effect_size,
                grad = contrast
            )
        }
    )
}


#' Standard-deviation of a set of coefficients as effect size
#'
#' This is intended as the effect size version of an ANOVA. For effect_sd, the effect size is the standard deviation of some coefficients about their mean. For effect_rss, it is the root sum of squared differences from the mean.
#'
#' @param coef The column numbers of the design matrix for the relevant coefficients.
#'
#' @return
#'
#' An object defining how to calculate an effect size.
#'
#' \code{effect_rss} may be better suited to comparing effect sizes from designs with differing numbers of coefficients, such as differential exon usage.
#'
#' @export
effect_sd <- function(coef) {
    n <- length(coef)
    assert_that(n > 1)   #n=2 case may be problematic

    list(
        signed = FALSE,
        df = n-1,

        calc = function(beta) sd(beta[coef]),

        constraint = function(effect_size) {
            target <- effect_size**2 * (length(coef)-1)
            function(beta) {
                bc <- beta[coef]
                mbc <- mean(bc)
                list(
                    score = target - sum((bc-mbc)**2),
                    grad = -2*(1-1/n)*(bc-mbc)
                )
            }
        }
    )
}


#' @rdname effect_sd
#' @export
effect_rss <- function(coef) {
    n <- length(coef)
    assert_that(n > 1)   #n=2 case may be problematic

    list(
        signed = FALSE,
        df = n-1,

        calc = function(beta) sqrt(sum( (beta[coef]-mean(beta[coef]))**2 )),

        constraint = function(effect_size) {
            target <- effect_size**2
            function(beta) {
                bc <- beta[coef]
                mbc <- mean(bc)
                list(
                    score = target - sum((bc-mbc)**2),
                    grad = -2*(1-1/n)*(bc-mbc)
                )
            }
        }
    )
}


#' Shift of mass effect
#'
#' Detect a "shift of mass" between two conditions. For example expression might move later in a time series in an experimental condition vs a control. If all expression shifted to a later time in the experimental condition, this would be given an effect size of 1. Conversely if all expression shifted to an earlier time, the effect size would be -1.
#'
#' @param coef1 Column numbers in the design matrix for the first condition, in some meaningful order.
#'
#' @param coef2 Corresponding column numbers for the second condition.
#'
#' \code{effect_shift_log2} is adapted to work with log2 scaled coefficients. This is almost certainly the version you want.
#'
#' @return
#'
#' An object defining how to calculate an effect size.
#'
#' @export
effect_shift <- function(coef1, coef2) {
    assert_that(length(coef1) == length(coef2))
    n <- length(coef1)
    sign_mat <- sign(outer(-seq_len(n),seq_len(n), `+`))

    list(
        signed = TRUE,
        df = 1,

        calc = function(beta) {
            beta1 <- beta[coef1]
            beta2 <- beta[coef2]
            a <- sum(outer(beta1,beta2)*sign_mat)
            b <- sum(beta1) * sum(beta2)
            a/b
        },

        constraint = function(effect_size) {
            function(beta) {
                beta1 <- beta[coef1]
                beta2 <- beta[coef2]
                a <- sum(outer(beta1,beta2)*sign_mat)
                b1 <- sum(beta1)
                b2 <- sum(beta2)
                grad <- rep(0.0, length(beta))
                grad[coef1] <-               colSums(beta2*t(sign_mat))/b1/b2 - a/b2/(b1*b1)
                grad[coef2] <- grad[coef2] + colSums(beta1*sign_mat)/b1/b2    - a/b1/(b2*b2)
                list(
                    score = a/b1/b2,
                    grad = grad
                )
            }
        }

        # constraint = function(effect_size) {
        #     weights <- sign_mat - effect_size
        #
        #     function(beta) {
        #         beta1 <- beta[coef1]
        #         beta2 <- beta[coef2]
        #         grad <- rep(0.0, length(beta))
        #         grad[coef1] <-               colSums(beta2*t(weights))
        #         grad[coef2] <- grad[coef2] + colSums(beta1*weights)
        #         list(
        #             score = sum(weights * outer(beta1,beta2)),
        #             grad = grad
        #         )
        #     }
        # }
    )
}


#' Adapt an effect size object to work with coefficients estimated from log2 transformed data or with a log2 link function
#'
#' Differential expression analysis usually fits coefficients to log transformed values, or with a log2 link function. This function adapts an effect size object to work with such log-scale coefficients.
#'
#' @param effect An object definining an effect size on coefficients with untransformed scale.
#'
#' @return
#'
#' An object defining how to calculate an effect size from coefficients on a logarithmic scale.
#'
#' @export
effect_link_log2 <- function(effect) {
    list(
        signed = effect$signed,
        df = effect$df,

        calc = function(beta) effect$calc(exp(beta*ln2)),

        constraint = function(effect_size) {
            cons <- effect$constraint(effect_size)

            function(beta) {
                e <- exp(beta*ln2)
                result <- cons(e)
                list(
                    score = result$score,
                    grad = ln2 * e * result$grad
                )
            }
        }
    )
}


#' @rdname effect_shift
#' @export
effect_shift_log2 <- function(coef1, coef2)
    effect_link_log2(effect_shift(coef1, coef2))


constrained_fit_newton <- function(y, X, devi, cons=NULL, offset=0, initial=NULL, dtol=1e-7, ctol=1e-7, btol=1e-7) {
    n <- length(y)
    m <- ncol(X)

    # Initial guess
    if (!is.null(initial))
        beta <- initial
    else {
        beta <- qr.coef(qr(X), devi$y_to_x(y) - offset)
        beta[is.na(beta)] <- 0.0
    }

    iters <- 0
    pred <- c(X %*% beta) + offset
    out <- devi$devi2(pred,y)
    value <- sum(out$devi)
    well_constrained <- FALSE
    repeat {
        C <- colSums(out$devi_x*X)
        Q <- t(X) %*% (out$devi_xx*X)

        last_well_constrained <- well_constrained
        well_constrained <- TRUE
        if (!is.null(cons)) {
            cons_out <- cons(beta)

            scale <- sqrt(sum(cons_out$grad*cons_out$grad))
            cstep <- cons_out$score / scale
            if (abs(cstep) > ctol) well_constrained <- FALSE

            C <- c(C, cstep)
            Q <- cbind(rbind(Q,cons_out$grad / scale),c(cons_out$grad / scale,0))
        }

        step <- solve(Q,C)[seq_len(m)]
        beta <- beta - step

        pred <- c(X %*% beta) + offset
        out <- devi$devi2(pred,y)
        old_value <- value
        value <- sum(out$devi)
        iters <- iters + 1

        if (last_well_constrained &&
            well_constrained &&
            all(abs(step) < btol) &&
            value <= old_value+dtol) break

        if (iters >= 100) {
            warning("Failed to converge after 100 iterations.")
            break
        }
    }

    list(beta=beta, deviance=value, iters=iters)
}



constrained_fit_cobyla <- function(y, X, devi, cons=NULL, offset=0, initial=NULL) {
    # Initial guess
    if (!is.null(initial))
        beta <- initial
    else {
        beta <- qr.coef(qr(X), devi$y_to_x(y) - offset)
        beta[is.na(beta)] <- 0.0
    }

    f <- function(beta)
        sum(devi$devi(c(X %*% beta) + offset, y))

    if (is.null(cons))
        g <- NULL
    else
        g <- function(beta) {
            s <- cons(beta)$score
            c(s, -s) #Equality constraint
        }

    result <- nloptr::cobyla(beta, fn=f, hin=g)

    if (result$convergence <= 0)
        warning("COBYLA failed to successfully converge.")

    list(beta=result$par, deviance=result$value, iters=result$iter)
}



constrained_fit_slsqp <- function(y, X, devi, cons=NULL, offset=0, initial=NULL) {
    # Initial guess
    if (!is.null(initial))
        beta <- initial
    else {
        beta <- qr.coef(qr(X), devi$y_to_x(y) - offset)
        beta[is.na(beta)] <- 0.0
    }

    f <- function(beta) {
        d <- devi$devi1(c(X %*% beta) + offset, y)
        list(objective=sum(d$devi), gradient=colSums(d$devi_x*X))
    }

    if (is.null(cons))
        g <- NULL
    else
        g <- function(beta) {
            cons_out <- cons(beta)
            list(constraints=cons_out$score, jacobian=cons_out$grad)
        }

    result <- nloptr::nloptr(beta, eval_f=f, eval_g_eq=g, opts=list(algorithm="NLOPT_LD_SLSQP", xtol_abs=1e-6))

    if (result$status <= 0)
        warning("SLSQP failed to successfully converge.")

    list(beta=result$solution, deviance=result$objective, iters=result$iter)
}


