# | -●- Topconfects 

(Under development) TOP results by CONfident efFECT Size. Topconfects is an R package intended for RNA-seq or microarray Differntial Expression analysis and similar, where we are interested in estimating many effect sizes -- one per gene -- from few samples.

Topconfects builds on [TREAT](http://bioinformatics.oxfordjournals.org/content/25/6/765.long) p-values offered by the limma and edgeR packages. It tries a range of fold changes, and uses this to rank genes by effect size while maintaining a given FDR. See [nest_confects](reference/nest_confects.html) for details.

<br/>

**Stop using p-values as a proxy for effect size.** The difference between a p-value of 1e-6 and 1e-9 has no practical meaning in terms of significance, however tiny p-values are often used as a proxy for effect size. This is a misuse, as they might simply reflect greater quality of evidence (for example RNA-seq average read count or microarray average spot intensity). It is better to reject a broader set of hypotheses, while maintaining sensible significance level. For example, a confidence interval rejects values of a parameter outside of that interval, at a given significance level. Here we will be concentrating on rejecting values of a parameter smaller than some size, while maintaining a False Discovery Rate.

**No need to guess the best fold change cutoff.** TREAT requires a fold change cutoff to be specified. The conclusion of the TREAT paper suggests manually adjusting this cutoff based on your data. Topconfects merely automates this: you specify a False Discovery Rate appropriate to your purpose, then read down the resulting ranked list of genes as far as you wish. The "confect" value given in the last row that you use is the fold change cutoff required for TREAT to produce that set of genes at the given FDR.

**Rank by an effect size that is meaningful.** Once we stop thinking in terms of zero vs non-zero effect and start thinking about effect size, the limitations of a linear or log-linear model become apparent. Topconfects therefore includes a variety of non-linear effect sizes.

* Examining the interaction of two factors, a log-linear model allows testing of odds ratios. However this will give a large effect size when, say, a proportion shifts from 0.01% to 0.1%. We may instead be interested in the difference of proportions, to avoid looking at such small shifts. Topconfects provides `effect_shift_log2` for this.

* In an ANOVA test, we need to test multiple coefficients, and we may want to treat each group equally rather than nominating a control group. Topconfects provides `effect_sd` for this, where the effect size is the standard deviation of a set of coefficients.

Topconfects contains its own implementation of `edgeR::glmTreat` which uses a constrained GLM fitting routine, in order to allow this. **This part is very much experimental!** The approximation that quasi-likelihood ratios are F-distributed seems quite robust, but that may be lost by adding non-linear constraints. A poorly chosen constraint might effectively shave off more degrees of freedom than intended.

<br/>

## Usage

Use [limma_confects](reference/limma_confects.html) or [edger_confects](reference/edger_confects.html) as part of your limma or edgeR analysis.

* [Example RNA-seq analysis](articles/fold_change.html)

* [Example looking for interaction using a non-(log)linear effect size](articles/nonlinear_effect.html)

<br/>

## Install

```r
install.packages("devtools")

devtools::install_github("pfh/topconfects")
```

<br/>

## Contact

Author: Paul Harrison [@paulfharrison](https://twitter.com/paulfharrison) paul.harrison@monash.edu

<br/>

## Future work

Gene-set enrichment tests. Here also the smallest p-value does not necessarily imply the greatest interest.

<br/>

## References

McCarthy, D. J., and Smyth, G. K. (2009). Testing significance relative to a fold-change threshold is a TREAT. *Bioinformatics* 25, 765-771. doi: 10.1093/bioinformatics/btp053 

