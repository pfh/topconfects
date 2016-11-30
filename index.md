# | (-● Topconfects 
(Under development) TOP results by CONfident efFECT Size. Topconfects is an R package intended for RNA-seq or microarray Differntial Expression analysis and similar, where we are interested in estimating many effect sizes -- one per gene -- from few samples.

Topconfects builds on [TREAT](http://bioinformatics.oxfordjournals.org/content/25/6/765.long) p-values offered by the limma and edgeR packages. It tries a range of fold changes, and uses this to rank genes by effect size while maintaining a given FDR. See [nest_confects](reference/nest_confects.html) for details.

<br/>

**Banish silly p-values and q-values.** The difference between a p-value of 1e-6 and 1e-9 has no practical meaning. What matters is effect sizes. If your p-values are this small, *you should be rejecting a broader set of hypotheses*. For example, a confidence interval rejects values of a parameter outside of that interval, at a given significance level. Here we will be concentrating on rejecting values of a parameter smaller than some size.

**No need to guess the fold change.** TREAT is awesome, but it needs a fold change cutoff to be specified, and this needs some prior knowledge of the biology involved. What you *can* specify is a False Discovery Rate appropriate to your purpose. Topconfects then gives a ranked list of the most interesting genes to look at by confident effect size.

<br/>

## Install

```r
install.packages("devtools")

devtools::install_github("pfh/topconfects")
```

## Usage

Use [limma_confects](reference/limma_confects.html) or [edger_confects](reference/edger_confects.html) as part of your limma or edgeR analysis.

* [Example](articles/fold_change.html)

## Contact

Author: Paul Harrison [@paulfharrison](https://twitter.com/paulfharrison) paul.harrison@monash.edu

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

