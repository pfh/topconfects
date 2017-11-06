# | -●- Topconfects 

TOP results by CONfident efFECT Size. Topconfects is an R package intended for RNA-seq or microarray Differntial Expression analysis and similar, where we are interested in placing confidence bounds on many effect sizes--one per gene--from few samples.

Topconfects builds on [TREAT](http://bioinformatics.oxfordjournals.org/content/25/6/765.long) p-values offered by the limma and edgeR packages. It tries a range of fold changes, and uses this to rank genes by effect size while maintaining a given FDR. This also produces confidence bounds on the fold changes, with adjustment for multiple testing. See [nest_confects](reference/nest_confects.html) for details.

<br/>

* **A principled way to avoid using p-values as a proxy for effect size.** The difference between a p-value of 1e-6 and 1e-9 has no practical meaning in terms of significance, however tiny p-values are often used as a proxy for effect size. This is a misuse, as they might simply reflect greater quality of evidence (for example RNA-seq average read count or microarray average spot intensity). It is better to reject a broader set of hypotheses, while maintaining a sensible significance level.

* **No need to guess the best fold change cutoff.** TREAT requires a fold change cutoff to be specified. Topconfects instead asks you specify a False Discovery Rate appropriate to your purpose. You can then read down the resulting ranked list of genes as far as you wish. The "confect" value given in the last row that you use is the fold change cutoff required for TREAT to produce that set of genes at the given FDR.

* **(experimental) Rank by a more meaningful effect size.** Once we stop thinking in terms of zero vs non-zero effect and start thinking about effect size, we have greater freedom to choose an effect size that is meaningful. For example, examining the interaction of two factors, a linear model on log expression levels allows testing of odds ratios. However this will give a large effect size when, say, a proportion shifts from 0.01% to 0.1%. We may instead be interested in the difference of proportions, to avoid looking at such small shifts. Topconfects provides [effect_shift_log2](reference/effect_shift.html) for this.

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

