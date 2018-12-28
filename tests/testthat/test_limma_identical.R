context("topconfects TREAT implementation should exactly match limma's")

library(limma)
library(topconfects)

set.seed(1234)
n <- 100
folds <- seq(-2,2,length.out=n)
row_means <- runif(n, min=0, max=5)
means <- outer(folds, c(-1,-1,1,1)) + row_means
sds <- rep(runif(n, min=0,max=2), 4)
mat <- rnorm(length(means), mean=means, sd=sds)
dim(mat) <- dim(means) 

design <- cbind(c(1,1,0,0), c(0,0,1,1))

fit <- lmFit(mat, design)
cfit <- contrasts.fit(fit, contrasts=c(-1,1))

# Using our version of limma's treat
confects_standard <- limma_confects(cfit, coef=1)
# Using limma's treat function
confects_limma <- topconfects:::limma_confects_limma(cfit, coef=1)

capture.output(file="test_limma_identical.out", {
    cat("\nUsing topconfects implementation\n")
    print(confects_standard)
    cat("\nUsing limma\n")
    print(confects_limma)
})

# Tables should match exactly
test_that("table identical", {
    expect_equal(confects_standard$table, confects_limma$table)
})

# p values should match exactly
test_that("p values identical", {
    for(lfc in c(0.5, 1.0, 2.0)) {
        p <- confects_standard$pfunc(1:n, lfc)
        p_limma <- confects_limma$pfunc(1:n, lfc)
        expect_equal(p, p_limma)
    }
})