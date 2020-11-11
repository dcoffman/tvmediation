# tvmediation <a><img src = 'man/figures/tvma_logo.png' align = "right" height = "125" /></a>
R package for fitting time-varying mediation models

## Overview

This package includes a set of functions for estimating mediation effects that vary over time. The package allows a time-varying mediator and time-varying continuous or binary outcome. The treatment variable is assumed to be time-invariant and may have either two or three levels. Confidence intervals for the indirect effect are obtained via bootstrap. The goal of this method is to assess whether the indirect effect varies as a function of time.

## Installation

To use the time varying mediation analysis package in R, you must first install the package and load it. Before that, make sure you have `R version 4.0.3` or greater. There are two ways to install the package from the CRAN (Comprehensive R Archive Network) repository, by using `install.packages` or the `devtools` function. 

```{r}
install.packages("tvmediation", dependencies = TRUE)
```
The equivalent code using `devtools` is:
```{r}
devtools::install_cran("tvmediation", dependencies = TRUE) 
# MAKE SURE YOU HAVE devtools INSTALLED
```
Alternatively, if you want to install the package directly from the github repository to access new or revised functions in development, the following code may be used:
```{r}
devtools::install_github("dcoffman/tvmediation", dependencies = TRUE) 
# MAKE SURE YOU HAVE devtools INSTALLED
```
## Getting started

```{r}
library(tvmediation)
```

## Getting help
Summarized versions of the function vignettes can be accessed through this [link](https://github.com/dcoffman/tvmediation/wiki).
If you encounter a bug, please file a minimal reproducible example on [github](https://github.com/dcoffman/tvmediation/issues).

