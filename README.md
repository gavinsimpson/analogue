# analogue

<!-- badges: start -->

[![R build
status](https://github.com/gavinsimpson/analogue/workflows/R-CMD-check/badge.svg)](https://github.com/gavinsimpson/analogue/actions)
[![codecov](https://codecov.io/gh/gavinsimpson/analogue/branch/master/graph/badge.svg?token=GG5NQfgRFu)](https://app.codecov.io/gh/gavinsimpson/analogue)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/analogue)](https://cran.r-project.org/package=analogue)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/analogue)](https://cran.r-project.org/package=analogue)
<!-- badges: end -->

## What is analogue?
**analogue** is an R package for use with palaeoecological data. Originally, **analogue** was intended as an R implementation of analogue methods such as analogue matching, <acronym title="Receiver Operator Characteristic">ROC</acronym> curves, and <acronym title="Modern Analogue Technique">MAT</acronym> transfer function models, and the computation of dissimilarity coefficients. Since then the scope of the package has grown to include a number of other methods applicable to data routinely encountered in palaeoecology and palaeolimnology.

### Features

 * Transfer functions
     * MAT
     * Weighted Averaging with monotonic, inverse, and classical deshrinking, with and without tolerance down-weighting
     * Principal Component Regression (using ecologically-relevant transformations)
     * Cross-validation (Bootstrapping, leave-one-out, *k*-fold)
     * Analogue statistics
 * Analogue matching
 * Dissimilarity coefficients
     * Chord, Bray-Curtis, Gower's Generalised coefficient, Manhattan, ...
 * Dissimilarity decisions thresholds
     * <acronym title="Receiver Operator Characteristic">ROC</acronym> curves
     * Monte Carlo resampling
     * Logistic regression
 * Stratigraphic diagrams
 * Principal curves

## Bugs, feature requests
Bug reports and feature requests should be filed as [issues](https://github.com/gavinsimpson/analogue/issues).

## Licence
analogue is released under the [GNU General Public Licence Version 2](https://www.r-project.org/Licenses/GPL-2).

## Links

 * [CRAN page](https://cran.r-project.org/package=analogue)
