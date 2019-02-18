# Topconfects 

TOP results by CONfident efFECT Size. Topconfects is intended for RNA-seq or microarray Differntial Expression analysis and similar, where we are interested in placing confidence bounds on many effect sizes--one per gene--from few samples.

## Topconfects is in Bioconductor 3.9

As at February 2019, Bioconductor 3.9 is the development version of Bioconductor. Future development will occur in the Bioconductor git server.

* [Topconfects in Bioconductor](https://bioconductor.org/packages/topconfects)


# Install for R 3.5

Once Bioconductor 3.9 is released (should be first half of 2019), please use the Bioconductor version, see above.

topconfects can be installed from the GitHub repository with:

```r
install.packages("devtools")

devtools::install_github("pfh/topconfects")

library(topconfects)
```


## Quasi-Likelihood non-linear effect size branch

This is an experimental branch. Rank by a more meaningful effect size: Once we stop thinking in terms of zero vs non-zero effect and start thinking about effect size, we have greater freedom to choose an effect size that is meaningful. For example, examining the interaction of two factors, a linear model on log expression levels allows testing of odds ratios. However this will give a large effect size when, say, a proportion shifts from 0.01% to 0.1%. We may instead be interested in the difference of proportions, to avoid looking at such small shifts.

This branch works, but is overly complicated for what it achieves. The current plan is to re-implement this idea using a simpler Wald test.

```r
devtools::install_github("pfh/topconfects@ql")

library(topconfectsql)
```


## Author

Topconfects is developed by Paul Harrison [@paulfharrison](https://twitter.com/paulfharrison) paul.harrison@monash.edu, for the [Monash Bioinformatics Platform](https://platforms.monash.edu/bioinformatics/).

This software is Copyright 2018 Paul Harrison, distributed under the terms of the GNU Lesser General Public License version 2.1 (see file LICENSE).

<a href="https://platforms.monash.edu/bioinformatics/"><img src="https://raw.githubusercontent.com/pfh/topconfects/master/MBP-logo.png" height="88"></a>


