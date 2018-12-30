
# Topconfects 

TOP results by CONfident efFECT Size. Topconfects is an R package intended for RNA-seq or microarray Differntial Expression analysis and similar, where we are interested in placing confidence bounds on many effect sizes---one per gene---from few samples, and ranking genes by these confident effect sizes.

Topconfects builds on [TREAT](http://bioinformatics.oxfordjournals.org/content/25/6/765.long) p-values offered by the limma and edgeR packages, or the "greaterAbs" test p-values offered by DESeq2. It tries a range of fold changes, and uses this to rank genes by effect size while maintaining a given FDR. This also produces confidence bounds on the fold changes, with adjustment for multiple testing.

* **A principled way to avoid using p-values as a proxy for effect size.** The difference between a p-value of 1e-6 and 1e-9 has no practical meaning in terms of significance, however tiny p-values are often used as a proxy for effect size. This is a misuse, as they might simply reflect greater quality of evidence (for example RNA-seq average read count or microarray average spot intensity). It is better to reject a broader set of hypotheses, while maintaining a sensible significance level.

* **No need to guess the best fold change cutoff.** TREAT requires a fold change cutoff to be specified. Topconfects instead asks you specify a False Discovery Rate appropriate to your purpose. You can then read down the resulting ranked list of genes as far as you wish. The "confect" value given in the last row that you use is the fold change cutoff required for TREAT to produce that set of genes at the given FDR.

The method is described in:

[Harrison PF, Pattison AD, Powell DR, Beilharz TH. 2018. Topconfects: a package for confident effect sizes in differential expression analysis provides improved usability ranking genes of interest. bioRxiv. doi:10.1101/343145](https://www.biorxiv.org/content/early/2018/06/11/343145)

<br/>

## Usage

Use [limma_confects](reference/limma_confects.html), [edger_confects](reference/edger_confects.html), or [deseq2_confects](reference/deseq2_confects.html) as part of your limma, edgeR, or DESeq2 analysis (the limma method is currently much faster than other methods).

* [Example RNA-seq analysis](articles/fold_change.html)

<br/>

## Install

Current git master branch specifies R 3.6. Please use release 1.1.4 if you are using R 3.5.

```r
install.packages("devtools")

devtools::install_github("pfh/topconfects@v1.1.4")

library(topconfects)
```

<br/>

## Quasi-Likelihood non-linear effect size branch

This is an experimental branch. Rank by a more meaningful effect size: Once we stop thinking in terms of zero vs non-zero effect and start thinking about effect size, we have greater freedom to choose an effect size that is meaningful. For example, examining the interaction of two factors, a linear model on log expression levels allows testing of odds ratios. However this will give a large effect size when, say, a proportion shifts from 0.01% to 0.1%. We may instead be interested in the difference of proportions, to avoid looking at such small shifts.

This branch works, but is overly complicated for what it achieves. The current plan is to re-implement this idea using a simpler Wald test.

```r
devtools::install_github("pfh/topconfects@ql")

library(topconfectsql)
```

* [Reference manual for topconfectsql](http://logarithmic.net/topconfects/ql/topconfectsql_1.0.1.pdf)
* [Example usage vignette](http://logarithmic.net/topconfects/ql/nonlinear_effect.html)

<br/>

## Contact

Topconfects is developed by Paul Harrison [@paulfharrison](https://twitter.com/paulfharrison) paul.harrison@monash.edu, for the [Monash Bioinformatics Platform](https://platforms.monash.edu/bioinformatics/).

<a href="https://platforms.monash.edu/bioinformatics/"><img src="https://raw.githubusercontent.com/pfh/topconfects/master/MBP-logo.png" height="88"></a>

<br/>

<!--
## Future work

Gene-set enrichment tests. Here also the smallest p-value does not necessarily imply the greatest interest.

<br/>
-->

## References

McCarthy, D. J., and Smyth, G. K. (2009). Testing significance relative to a fold-change threshold is a TREAT. *Bioinformatics* 25, 765-771. doi: 10.1093/bioinformatics/btp053 

