# analogue

[![Build Status](https://travis-ci.org/gavinsimpson/analogue.svg?branch=master)](https://travis-ci.org/gavinsimpson/analogue)  [![Build status](https://ci.appveyor.com/api/projects/status/hc8dbxrim2nj3c1i/branch/master)](https://ci.appveyor.com/project/gavinsimpson/analogue/branch/master) [![CRAN version](http://www.r-pkg.org/badges/version/analogue)](http://cran.rstudio.com/web/packages/analogue/index.html) [![](http://analogue.r-pkg.org/badges/grand-total/analogue)](http://cran.rstudio.com/web/packages/analogue/index.html)

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
analogue is released under the [GNU General Public Licence Version 2](http://www.gnu.org/licenses/gpl-2.0.html).

## Links

 * [CRAN page](http://cran.r-project.org/package=analogue)
