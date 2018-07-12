
# Nonlinear group effect size definitions
#
# Once supplied with parameters and column names, these take a matrix [feature,coef] of coefficients and produce a list(effect=, gradient=). effect is a single number. gradient is a matrix [feature,coef].
#


group_effect_unlog2 <- 
function(group_effect) 
function(col_names) {
    group_effect <- group_effect(col_names)

    list(
        limits=group_effect$limits,

        calc=function(beta) {
            exp2_beta <- exp(beta*ln2)
            inner <- group_effect$calc(exp2_beta)
            list(
                effect = inner$effect,
                gradient = ln2*exp2_beta*inner$gradient)
        }
    )
}

group_effect_shift <-
function(coef1, coef2)
function(col_names) {
    coef1 <- resolve_coef(coef1, col_names)
    coef2 <- resolve_coef(coef2, col_names)
    assert_that(length(coef1) == 1)
    assert_that(length(coef2) == 1)

    n <- length(col_names)

    list(
        limits=c(-1,1),

        calc=function(beta) {
            m <- nrow(beta)
            assert_that(m >= 2, msg="groups must contain at least 2 members")

            num <- 0.0
            den <- 0.0
            grad_num <- matrix(0, nrow=m,ncol=n)
            grad_den <- matrix(0, nrow=m,ncol=n)
            for(i in seq_len(m)) {
                for(j in seq_len(m)) {
                    s <- sign(j-i)
                    mul <- beta[i,coef1]*beta[j,coef2]
                    num <- num + s*mul
                    grad_num[i,coef1] <- grad_num[i,coef1] + s*beta[j,coef2]
                    grad_num[j,coef2] <- grad_num[j,coef2] + s*beta[i,coef1]
                    den <- den + mul
                    grad_den[i,coef1] <- grad_den[i,coef1] + beta[j,coef2]
                    grad_den[j,coef2] <- grad_den[j,coef2] + beta[i,coef1]
                }
            }

            list(
                effect=num/den,
                gradient=(grad_num*den-grad_den*num)/ (den*den)
            )
        }
    )
}

group_effect_shift_unlog2 <- function(coef1, coef2)
    group_effect_unlog2(group_effect_shift(coef1, coef2))


