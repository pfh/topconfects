
#' @export
edger_confects <- function(fit, coef=ncol(fit$design), contrast=NULL, fdr=0.05, max=30.0, step=0.05) {
    assert_that(is(fit, "DGEGLM"))

    n <- nrow(fit)

    pfunc <- function(i, mag) {
        top_treats <-
            glmTreat(fit, coef=coef, contrast=contrast, lfc=mag) %>%
            topTags(n=n, sort="none")
        top_treats$table$PValue[i]
    }

    confects <- nest_confects(n, pfunc, fdr=fdr, max=max, step=step)

    top_tags <- glmQLFTest(efit,coef=coef,contrast=contrast) %>% topTags(n=n,sort="none")
    logFC <- top_tags$table$logFC[confects$index]
    confects$signed_confect <- sign(logFC) * confects$confect
    confects$logFC <- logFC
    confects$logCPM <- top_tags$table$logCPM[confects$index]
    confects$name <- rownames(top_tags$table)[confects$index]

    confects
}
