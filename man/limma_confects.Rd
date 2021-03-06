% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/limma.R
\name{limma_confects}
\alias{limma_confects}
\title{Confident log2 fold changes based on a limma fit object}
\usage{
limma_confects(
  fit,
  coef = NULL,
  fdr = 0.05,
  step = 0.001,
  trend = FALSE,
  full = FALSE
)
}
\arguments{
\item{fit}{A limma MArrayLM object.}

\item{coef}{Number or name of coefficient or contrast to test.}

\item{fdr}{False Discovery Rate to control for.}

\item{step}{Granularity of log2 fold changes to test.}

\item{trend}{Should \code{eBayes(trend=TRUE)} be used?}

\item{full}{Include some further statistics used to calculate confects in the
output, and also include FDR-adjusted p-values that effect size is non-zero
(note that this is against the spirit of the topconfects approach).}
}
\value{
See \code{\link{nest_confects}} for details of how to interpret the result.
}
\description{
For all possible absolute log2 fold changes (LFC), which genes have at least
this fold change at a specified False Discovery Rate (FDR)?
}
\details{
Results are presented in a table such that for any given LFC, if the reader
chooses the genes with \code{abs(confect)} less than this they are assured
that this set of genes has at least this LFC (with the specified FDR). Once
this set of genes is selected, the confect values provide confidence bounds
with False Coverage-statement Rate at the same level as the FDR.

\code{fit} should be produced using \code{lmFit}. It is not necessary to use
\code{eBayes}, this function calls \code{eBayes} itself.

To test contrasts, this function can also be used with the result of
\code{contrasts.fit}, but limma's handling of weights may be approximate (for
example if \code{voom} has been used). For exact results for a contrast, use
\code{contrastToCoef} to adjust the design matrix given to \code{lmFit}.
}
\examples{

#Prepare a data set
library(NBPSeq)
library(edgeR)
library(limma)
data(arab)
dgelist <- DGEList(arab)
dgelist <- calcNormFactors(dgelist)
cpms <- cpm(dgelist, log=TRUE)
# Retain genes with more than a geometric mean of 2 RPM
# (about 5 reads per sample)
cpms <- cpms[rowMeans(cpms) >= 1,]

# Fit linear model for each gene
treatment <- c(FALSE,FALSE,FALSE,TRUE,TRUE,TRUE)
batch <- factor(c(1,2,3,1,2,3))
design <- model.matrix(~ treatment + batch)
fit <- lmFit(cpms, design)

# Calculate top confects
# As voom has not been used, it is necessary to use trend=TRUE
limma_confects(fit, "treatmentTRUE", trend=TRUE)

}
