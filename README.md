# VSALM: Variance-Smoothing Area Level Models for Small Area Estimation of Indicators

## Overview 

This package contains Stan code and R code to fit area level models for small area estimation of proportions.

## Install

`VSALM` uses Stan to conduct Bayesian inference. Before installing `VSALM`, install `rstan`, making sure that a functioning C++ compiler is installed as well. Once this is complete, `VSALM` may be installed from GitHub.


```
install.packages("devtools")
devtools::install_github("peteragao/VSALM", dependencies=TRUE)
```

## To do

1. Add tests
2. Add vignette
