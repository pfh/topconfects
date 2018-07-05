
# Non-linear effect size definitions

ln2 <- log(2)



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


#' Adapt an effect object to work with coefficients estimated from log2 transformed data or with a log2 link function
#'
#' Differential expression analysis usually fits coefficients to log2 transformed values, or with a log2 link function. This function adapts an effect size object to work with such log2-scale coefficients.
#'
#' @param effect An object definining an effect size on coefficients with untransformed scale.
#'
#' @return
#'
#' An object defining how to calculate an effect size from coefficients on a logarithmic scale.
#'
#' @export
effect_link_log2 <- 
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
#' @exporteffect_contrast <- 
function(contrast) 
function(col_names) {
    assert_that(is.numeric(contrast))

    if (!is.null(names(effect))) {
        assert_that(all( names(contrast) %in% col_names ))

        new_contrast <- rep(0.0, ncol(design))
        new_contrast[ match(names(contrast), col_names) ] <- contrast
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


