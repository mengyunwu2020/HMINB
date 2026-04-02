# HMINB

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**HMINB** is an R package for fitting a Heterogeneous Multiple-Inflated Negative Binomial model to count data with possible inflation at multiple points. The model automatically selects inflation points and performs variable selection in both the count and inflation components.

## Features

- Handles multiple inflated values with individual-specific inflation probabilities via cumulative logit.
- Fused LASSO penalty to automatically identify true inflation points.
- Adaptive LASSO for variable selection in NB and inflation components.
- EM algorithm with efficient updates.

## Installation

Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("mengyunwu2020/HMINB")
