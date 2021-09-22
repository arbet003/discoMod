
<!-- README.md is generated from README.Rmd. Please edit that file -->

# discoMod

<!-- badges: start -->

<!-- badges: end -->

The goal of discoMod is to identify “modules” (subgroups of correlated
genes or other ’omics features), and then test whether the network
structure within a module differs between two phenotype groups. By
default, the network struture is measured using Spearman’s correlation,
but other options are provided. Here is a brief overview of the three
main functions:

  - *find\_modules()*: model-based clustering to identify modules of
    correlated genes
  - *test\_modules()*: test whether the network structure within a
    module differs between the two phenotype groups
  - *corrheatmap()*: for a given module, plot correlation heatmaps to
    visualize differences between the two phenotype groups

## Installation

Use the following code to install the discoMod R package:

``` r
# install.packages("devtools")
devtools::install_github("arbet003/discoMod")
```

## Tutorial

For a tutorial, see the example near the end of the following help page:

``` r
library(discoMod)
?test_modules
```

## Reference

When citing the discoMod R package, please use the following:

Arbet, J., Zhuang, Y., Litkowski, E., Saba, L., & Kechris, K. (2021). Comparing Statistical Tests for Differential Network Analysis of Gene Modules. Frontiers in genetics, 12, 748.
