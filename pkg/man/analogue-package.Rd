\name{analogue-package}
\alias{analogue-package}
\alias{analogue}
\docType{package}
\title{
Analogue and weighted averaging methods for palaeoecology
}
\description{
Fits Modern Analogue Technique and Weighted Averaging transfer 
  	     function models for prediction of environmental data from species 
	     data, and related methods used in palaeoecology.
}
\details{
\tabular{ll}{
Package: \tab analogue\cr
Type: \tab Package\cr
Version: \tab 0.11-99\cr
Date: \tab $Date: 2013-10-05 15:11:52 -0600 (Sat, 05 Oct 2013) $\cr
Depends: \tab R (>= 2.15.0), vegan (>= 1.17-12), lattice, rgl\cr
Imports: \tab mgcv, MASS, stats, graphics, grid, brglm, princurve\cr
Suggests: \tab testthat\cr
Authors@R: \tab c(person(given = c("Gavin", "L."), family = "Simpson",
                    role = c("aut", "cre"),
                    email = "ucfagls@gmail.com"),
             person(given = "Jari", family = "Oksanen",
                    role = "aut"))\cr
BugReports: \tab 
http://r-forge.r-project.org/tracker/?func=browse&group_id=69&atid=338\cr
NeedsCompilation: \tab yes\cr
License: \tab GPL-2\cr
ByteCompile: \tab true\cr
URL: \tab http://analogue.r-forge.r-project.org\cr
Packaged: \tab 2013-12-13 05:23:57 UTC; gavin\cr
Built: \tab R 3.0.2; x86_64-unknown-linux-gnu; 2013-12-13 05:27:27 UTC; unix\cr
}

Index:
\preformatted{
ImbrieKipp              Imbrie and Kipp foraminifera training set
Pollen                  North American Modern Pollen Database
RMSEP                   Root mean square error of prediction
Stratiplot              Palaeoecological stratigraphic diagrams
abernethy               Abernethy Forest Pollen Sequence
analog                  Analogue matching
analogue-package        Analogue and weighted averaging methods for
                        palaeoecology
bayesF                  Bayes factors
bootstrap               Bootstrap estimation and errors
bootstrap.wa            Bootstrap estimation and errors for WA models
bootstrapObject         Bootstrap object description
caterpillarPlot         Caterpillar plot of species' WA optima and
                        tolerance range.
chooseTaxa              Select taxa (variables) on basis of maximum
                        abundance attained and number of occurrences.
cma                     Close modern analogues
crossval                Cross-validation of palaeoecological transfer
                        function models
densityplot.residLen    Lattice density plot for residual lengths
deshrink                Deshrinking techniques for WA transfer
                        functions
dissimilarities         Extract dissimilarity coefficients from models
distance                Flexibly calculate dissimilarity or distance
                        measures
fitted.logitreg         Fitted values for the training set from
                        logistic regression models
fuse                    Fused dissimilarities
getK                    Extract and set the number of analogues
gradientDist            Positions of samples along a unit-length
                        ordination gradient.
hist.residLen           Histogram plot for residual lengths
histogram.residLen      Lattice histogram plot for residual lengths
join                    Merge species data sets on common columns
                        (species)
logitreg                Logistic regression models for assessing
                        analogues/non-analogues
mat                     Modern Analogue Technique transfer function
                        models
mcarlo                  Monte Carlo simulation of dissimilarities
minDC                   Extract minimum dissimilarities
optima                  Weighted averaging optima and tolerance ranges
panel.Loess             Loess smooths to stratigraphic diagrams
panel.Stratiplot        Panel function for stratigraphic diagrams
pcr                     Prinicpal component regression transfer
                        function models
performance             Transfer function model performance statistics
plot.dissimilarities    Plots the distribution of extracted
                        dissimilarities
plot.logitreg           Produces plots of analogue logistic regression
                        models
plot.mat                Plot diagnostics for a mat object
plot.mcarlo             Plot Monte Carlo simulated dissimilarity
                        distributions
plot.minDC              Plot of minimum dissimilarity per sample
plot.prcurve            Plot a fitted principal curve in PCA space
plot.residLen           Plot method for residual lengths
plot.roc                Plot ROC curves and associated diagnostics
plot.wa                 Plot diagnostics for a weighted averaging model
plot3d.prcurve          Interactive 3D plof of a principal curve in
                        principal coordinate space
prcurve                 Fits a principal curve to m-dimensional data
predict.logitreg        Posterior probability of analogue-ness for
                        fossil samples
predict.mat             Predict method for Modern Analogue Technique
                        models
predict.pcr             Predicted values from a principal components
                        regression
predict.wa              Predict from a weighted average model
rankDC                  Rank correlation between environmental and
                        species dissimilarities.
reconPlot               Stratigraphic plots of palaeoenvironmental
                        reconstructions
residLen                Squared residual length diagnostics
residuals.prcurve       Residuals of a principal curve fit.
rlgh                    Round Loch of Glenhead Diatoms
roc                     ROC curve analysis
scores.prcurve          'scores' method for principal curve objects of
                        class '"prcurve"'.
screeplot.mat           Screeplots of model results
smoothSpline            Smoother plugin function for use in fitting a
                        principal curve
splitSample             Select samples from along an environmental
                        gradient
sppResponse             Species responses along gradients.
stdError                Standard error of MAT fitted and predicted
                        values
summary.analog          Summarise analogue matching results
summary.bootstrap.mat   Summarise bootstrap resampling for MAT models
summary.cma             Summarise the extraction of close modern
                        analogues
summary.mat             Summarise Modern Analogue Technique models
summary.predict.mat     Summarise MAT model predictions
swapdiat                SWAP sub-fossil diatom and pH training set
swappH                  SWAP sub-fossil diatom and pH training set
timetrack               Timetracks of change in species composition
tran                    Common data transformations and
                        standardizations
varExpl                 Variance explained by ordination axes
wa                      Weighted averaging transfer functions
weightedCor             Weighted correlation test of WA reconstruction
}

Further information is available in the following vignettes:
\tabular{ll}{
\code{analogue_methods} \tab Analogue Methods in Palaeoecology (source, pdf)\cr
}
}
\author{
Gavin L. Simpson, Jari Oksanen

Maintainer: Gavin L. Simpson <ucfagls@gmail.com>
}
\keyword{ package }
