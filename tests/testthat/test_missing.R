context("Test missing values handled correctly")

library(topconfects)


effect <- c(1, NA, 10, 3, -20, NA, 7, NA)
se     <- c(1, 0,  2,  2,  1,  0, NA, 1)
good   <- !is.na(effect) & !is.na(se)


res1 <- normal_confects(effect, se, full=TRUE)
sub <- res1$table[ good[res1$table$index], ]

res2 <- normal_confects(effect[good], se[good], full=TRUE)

test_that("results_match", {
    expect_equal(sub$confect, res2$table$confect)
    expect_equal(sub$effect, res2$table$effect)
    expect_equal(sub$se, res2$table$se)
})

