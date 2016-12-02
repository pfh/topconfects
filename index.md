# | (-● Topconfects 

(Under development) TOP results by CONfident efFECT Size. Topconfects is an R package intended for RNA-seq or microarray Differntial Expression analysis and similar, where we are interested in estimating many effect sizes -- one per gene -- from few samples.

Topconfects builds on [TREAT](http://bioinformatics.oxfordjournals.org/content/25/6/765.long) p-values offered by the limma and edgeR packages. It tries a range of fold changes, and uses this to rank genes by effect size while maintaining a given FDR. See [nest_confects](reference/nest_confects.html) for details.

<br/>

**Stop using p-values as a proxy for effect size.** The difference between a p-value of 1e-6 and 1e-9 has no practical meaning in terms of significance, however tiny p-values are often used as a proxy for effect size. This is a misuse, as they might simply reflect a greater amount of evidence (for example RNA-seq average read count or microarray average spot intensity). It is better to reject a broader set of hypotheses, while maintaining sensible significance level. For example, a confidence interval rejects values of a parameter outside of that interval, at a given significance level. Here we will be concentrating on rejecting values of a parameter smaller than some size, while maintaining a False Discovery Rate.

**No need to guess the best fold change cutoff.** TREAT requires a fold change cutoff to be specified. The conclusion of the TREAT paper suggests manually adjusting this cutoff based on your data. Topconfects merely automates this: you specify a False Discovery Rate appropriate to your purpose, then read down the resulting ranked list of genes as far as you wish. The "confect" value given in the last row that you use is the fold change cutoff required for TREAT to produce that set of genes at the given FDR.

<br/>

## Usage

Use [limma_confects](reference/limma_confects.html) or [edger_confects](reference/edger_confects.html) as part of your limma or edgeR analysis.

* [Example RNA-seq analysis](articles/fold_change.html)

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

```
           Factor 2
             -  +
           .------  
          -| a  b
Factor 1   |
          +| c  d
```

In a two-factor experiment where the interaction of factors is of interest, the interaction term of a linear model may not be a good way to quantify the effect size. In Differential Expression analysis, this will be a log odds ratio, log((b/a)/(d/c)), however a more interesting value might be the difference in proportions, b/(a+b)-d/(c+d).

Differential exon usage is an instance of this, in which exon becomes one of the factors in the experiment.

Another line for future work is effect sizes for gene-set enrichment tests. The smallest p-value does not necessarily imply the greatest interest.

<br/>

## References

McCarthy, D. J., and Smyth, G. K. (2009). Testing significance relative to a fold-change threshold is a TREAT. *Bioinformatics* 25, 765-771. doi: 10.1093/bioinformatics/btp053 

