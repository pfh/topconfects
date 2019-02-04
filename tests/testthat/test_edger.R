context("Test edgeR-based routines")

# test_edger.out can be examined

library(edgeR)
library(topconfects)

set.seed(1234)

n <- 100
folds <- seq(-4,4,length.out=n)
row_means <- runif(n, min=0, max=5)
lib_scale <- c(1,2,3,4)
means <- 2^(outer(folds, c(-0.5,-0.5,0.5,0.5))) * row_means * rep(lib_scale,each=n)
counts <- rnbinom(length(means), mu=means, size=1/0.1)
dim(counts) <- dim(means) 

design <- cbind(c(1,1,0,0), c(0,0,1,1))

y <- DGEList(counts)
y <- calcNormFactors(y)
y <- estimateDisp(y, design)

fit <- glmQLFit(y, design)

set.seed(12345)
confects <- edger_confects(fit, contrast=c(-1,1))

achieved_fcrs <- cumsum(folds[confects$table$index]/confects$table$confect <= 1) / (1:n)

capture.output(file="test_edger.out", {
    cat("\ncounts\n")
    print(counts)
    cat("\ncolSums(counts)\n")
    print(colSums(counts))
    cat("\ndf.prior\n")
    print(fit$df.prior)
    cat("\nconfects\n")
    print(confects)
    cat("\nconfects$table\n")
    print(confects$table)
    cat("\nAchieved FCRs\n")
    print(achieved_fcrs)
})

test_that("FCR of 5% for all heads of table", {
    expect(all(achieved_fcrs <= 0.05, na.rm=TRUE), "FCR not controlled")
})
