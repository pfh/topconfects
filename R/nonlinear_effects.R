
# Non-linear effect size definitions

ln2 <- log(2)

resolve_coef <- function(coef, col_names) {
    if (is.character(coef))
        coef <- match(coef, col_names)

    assert_that(is.numeric(coef))
    assert_that(!any(is.na(coef)))
    assert_that(all(coef %in% seq_along(col_names)))

    coef
}


#' @export
effect_expr <- 
function(formula, df=1, limits=NULL)
function(col_names) {
    f <- deriv(formula, namevec=col_names)

    list(
        df=df,
        limits=limits,

        calc = function(beta) {
            env <- as.list(beta)
            names(env) <- col_names
            result <- eval(f, env)
            list(
                effect = as.vector(result),
                gradient = as.vector(attr(result,"gradient")))
        }
    )
}


#' Adapt an effect object to work with coefficients estimated from log2 transformed data
#'
#' Differential expression analysis usually fits coefficients to log2 transformed values. This function adapts an effect size object to work with such log2-scale coefficients.
#'
#' @param effect An object definining an effect size on coefficients with untransformed scale.
#'
#' @return
#'
#' An object defining how to calculate an effect size from coefficients on a logarithmic scale.
#'
#' @export
effect_unlog2 <- 
function(effect) 
function(col_names) {
    effect <- effect(col_names)
    list(
        df = effect$df,
        limits = effect$limits,

        calc = function(beta) {
            exp2_beta <- exp(beta*ln2)
            inner <- effect$calc(exp2_beta)
            list(
                effect = inner$effect,
                gradient = ln2*exp2_beta*inner$gradient)
        }
    )
}



#' Linear contrast effect.
#'
#' This is included for completeness and testing. limma_confects() provides a faster version.
#'
#' @param contrast Either a named vector with names matching the design matrix, or a vector of the same length as the number of columns in the design matrix. The effect size is the dot product of this with the fitted coefficients.
#'
#' @return
#'
#' An object defining how to calculate an effect size.
#'
#' @export
effect_contrast <- 
function(contrast) 
function(col_names) {
    assert_that(is.numeric(contrast))

    if (!is.null(names(contrast))) {
        new_contrast <- rep(0.0, ncol(design))
        new_contrast[ resolve_coef(names(contrast), col_names) ] <- contrast
        contrast <- new_contrast
    }

    list(
        df = 1,

        calc = function(beta) {
            list(
                effect = sum(beta * contrast),
                gradient = contrast)
        }
    )
}


#' @export
effect_shift <-
function(coef1, coef2)
function(col_names) {
    coef1 <- resolve_coef(coef1, col_names)
    coef2 <- resolve_coef(coef2, col_names)
    assert_that(length(coef1) == length(coef2))

    n <- length(col_names)
    m <- length(coef1)
    assert_that(m > 1)

    list(
        df = 1,
        limits = c(-1,1),

        calc = function(beta) {
            num <- 0.0
            den <- 0.0
            grad_num <- rep(0.0, n)
            grad_den <- rep(0.0, n)
            for(i in seq_len(m)) {
                for(j in seq_len(m)) {
                    s <- sign(j-i)
                    mul <- beta[coef1[i]]*beta[coef2[j]]
                    num <- num + s*mul
                    grad_num[coef1[i]] <- grad_num[coef1[i]] + s*beta[coef2[j]]
                    grad_num[coef2[j]] <- grad_num[coef2[j]] + s*beta[coef1[i]]
                    den <- den + mul
                    grad_den[coef1[i]] <- grad_den[coef1[i]] + beta[coef2[j]]
                    grad_den[coef2[j]] <- grad_den[coef2[j]] + beta[coef1[i]]
                }
            }

            list(
                effect=num/den,
                gradient=(grad_num*den-grad_den*num)/ (den*den)
            )
        }
    )
}


#' @export
effect_shift_unlog2 <- function(coef1, coef2)
    effect_unlog2(effect_shift(coef1, coef2))




