
<!-- README.md is generated from README.Rmd. Please edit that file -->

# distrom

<!-- badges: start -->

[![CRAN Version](https://www.r-pkg.org/badges/version/distrom)](https://www.r-pkg.org/pkg/distrom)
[![CRAN Posit Mirror Downloads](https://cranlogs.r-pkg.org/badges/grand-total/distrom)](https://www.r-pkg.org/pkg/distrom)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-success.svg)](https://lifecycle.r-lib.org/articles/stages.html)

<!-- badges: end -->

The R package *distrom* contains functions for computing a distributed
multinomial regression. The main function is `dmr()` which takes a
matrix of `covars` and a matrix of multinomial `counts` as input.
Independent Poisson log regressions of the form `counts ~ covars` are
then fit for each multinomial count. These independent Poisson log
regressions are estimated in parallel using the `parallel` and `gamlr`
packages which allows for easy in-memory parallelization and
distribution across multiple machines. This parallelization is essential
for use cases such as text analysis where the `counts` matrix consists
of many tokenized documents and can grow to billions of observations. In
the text analysis use case token counts are modeled as arising from a
multinomial distribution that is dependent upon the article attributes
contained in the `covars` matrix.

To cite this package, use “Taddy (2015), Distributed Multinomial
Regression, Annals of Applied Statistics”.

## Links

For a description of the functions in the *distrom* package please read
the reference manual: [distrom
manual](https://cran.r-project.org/web/packages/distrom/distrom.pdf)

For a detailed explanation of distributed multinomial regression and
example use cases see: [*Taddy (2015), Distributed Multinomial
Regression, Annals of Applied
Statistics*](https://arxiv.org/abs/1311.6139)

For information on the related *gamlr* package please read the [gamlr
manual](https://cran.r-project.org/web/packages/gamlr/gamlr.pdf) or
visit the [gamlr repository](https://github.com/TaddyLab/gamlr).

For information on the related *textir* package please read the [textir
manual](https://cran.r-project.org/web/packages/textir/textir.pdf) or
visit the [textir repository](https://github.com/TaddyLab/textir).

## Installation

To install the **stable** version from
[CRAN](https://cran.r-project.org/package=distrom):

``` r
install.packages("distrom")
```

To install the **development** version from
[GitHub](https://github.com/TaddyLab/distrom):

``` r
# install.packages("remotes")
remotes::install_github("TaddyLab/distrom")
```
