# tvmediation <a><img src = 'man/figures/tvma_logo.png' align = "right" height = "125" /></a>
R package for fitting time-varying mediation models

## Overview

This package includes a set of functions to measure the effect of mediating variables over time. The goal of this method is to take into consideration the dynamic effect that a mediator may have over time; thus the realistic nature of an indirect effect (if any) over time is characterized as opposed to the traditional approach of estimating an effect based on its single occurrence.

## Installation

To use the time varying mediation analysis package in R, you must first install the package and load it. Before that, make sure you have `R version 4.0.2`. There are two ways to install the package from the CRAN (Comprehensive R Archive Network) repository, by using `install.packages` or `devtools` function. 

```{r}
install.packages("tvmediation", dependencies = TRUE)
```
The equivalent code using `devtools` is:
```{r}
devtools::install_cran("tvmediation", dependencies = TRUE) # MAKE SURE YOU HAVE devtools INSTALLED AND LOADED #
```
Alternatively, if you want to install the package directly from the github repository to access new or revised functions in development, the following code is required:
```{r}
devtools::install_github("dcoffman/tvmediation", dependencies = TRUE) # MAKE SURE YOU HAVE devtools INSTALLED AND LOADED #
```
## Getting started

```{r}
library(tvmediation)
```

## Getting help
Summarized versions of the function vignettes can be accessed from this [link](https://github.com/dcoffman/tvmediation/wiki).
If you encounter a clear bug, please file a minimal reproducible example on [github](https://github.com/dcoffman/tvmediation/issues).

