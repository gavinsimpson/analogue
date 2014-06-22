\name{NEWS}
\title{Latest changes to the \pkg{analogue} package}
\encoding{UTF-8}

\section{Changes in version 0.14-0}{

  \subsection{General summary}{

    \itemize{
      \item Version 0.14-0 is less of an update than that provided by
      version 0.12-0 and sees a continuation of the changes to
      dependencies to have packaes in Imports rather than Depends.

      \item The main development of \pkg{analogue} now takes place on
      \href{https://github.com/gavinsimpson/analogue/}{github} and bugs
      and feature requests should be posted there. The Travis continuous
      integration system is used to automatically check the package as
      new code is checked in.

      \item There are several new functions and methods and a few bug
      fixes, the details of which are given below.

      \item The main user-visible change over 0.12-0 is the deprecation
      of the \code{plot3d.prcurve} method. The functionality is now in
      new function \code{Plot3d} and \code{plot3d.prcurve} is deprecated
      and if called needs to use the full function name. This change is
      to make \pkg{analogue} easier to install on MacOS X as now
      \pkg{rgl} is not needed to install analogue. If you want to plot
      the principal curve in an interactive 3d view, you'll need to get
      \pkg{rgl} installed first.
    }

  } % general

  \subsection{New Features}{

    \itemize{
      \item \code{n2} is a new utility function to calculate Hill's N2
      for sites and species.

      \item \code{optima} can now compute bootstrap WA optima and
      uncertainty.

      \item \code{performance} has a new method for objects of class
      \code{crossval}.

      \item \code{timetrack} had several improvements including a new
      \code{precict} method, which allows further points to be added to
      an existing timetrack, a \code{points} method to allow the
      addition of data to an existing timetrack plot, and the
      \code{plot} method can create a blank plotting region allowing
      greater customisation.

      \item \code{prcurve} gets \code{predict} and \code{fitted} methods
      to predict locations of new samples on the principal curve and
      extract the locations of the training samples respectively.
    }

  } % new features

  \subsection{Bug Fixes}{

    \itemize{
      \item The calculation of AUC in \code{roc} wasn't working
      correctly in some circumstances with just a couple of groups.

      \item \code{crossval.pcr} had a number of bugs in the k-fold CV
      routine which were leading to errors and the function not working.

      The progress bar was not being updated correctly either.

      \item \code{predict.pcr} was setting \code{ncomp} incorrectly if
      not supplied by the user.

      \item \code{ChiSquare} wasn't returning the transformation
      parameters required to transform leave-out data during
      crossvalidation or new samples for which predictions were
      required.

      \item \code{plot3d.prcurve} was not using the \code{data} and
      \code{ordination} components of the returned object. Note this
      function is now deprecated; see earlier note above.

      \item \code{predict.pcr} was incorrectly calling the internal
      function \code{fitPCR} with the \code{:::} operator.
  
    }
    
  }

  \subsection{DEPRECATED}{
    
    \itemize{
      \item \code{plot3d.prcurve} is deprecated. Functionality is in new
      function \code{Plot3d}.
    }
    
  } % deprecated
}
  
\section{Changes in version 0.12-0}{

  \subsection{General summary}{

    \itemize{
      \item Version 0.12-0 represents a major update of \pkg{analogue},
      including the addition of new functionality and long-required
      improvements to the computation of dissimilarity matrices which
      now uses faster C code.

      \item Several of the dependencies are now only imported into the
      namespace for \pkg{analogue}. As a result, existing code that did
      not already load these packages will need to be updated to load
      them via a call to \code{library()}.

      \item Highlights of the release are listed below. For a full
      record of changes made during the development leading up to this
      release see the ChangeLog in the package sources.
    }

  } % general

  \subsection{New Features}{

    \itemize{
      \item Dissimilarities are now computed using fast C code in
      function \code{distance()}. The new functions are tested against
      the old code (now available as function \code{oldDistance()}) to
      ensure the new code performs correctly.

      \item Fitting principal curves to palaeo data and working with the
      fitted objects is now much easier.

      \code{prcurve()} now returns each of the fitted smooth models as
      part of component \code{smooths}, allowing further examination of
      the the individual fits.

      There are now methods for \code{lines()} and \code{points} to
      allow greater flexibility in producing plots of the fitted
      principal curve. A \code{scores()} method to extract ordination
      scores for prcinipal curves is also available. A \code{resid()}
      method is also available to return various types of residuals for
      the principal curve.

      \code{smoothGAM()} is a new plugin smoother function for use with
      \code{prcurve()}, which allows GAMs to be used as smoothers fitted
      to individual variables. It is significantly slower than
      \code{smoothSpline()}, but can handle non-Gaussian responses. The
      underlying models are fitted using \code{gam()} from \pkg{mgcv}.

      As principal curves are smooth curves in high-dimensions it makes
      sense to view them in 3D. \pkg{analogue} has a function for that:
      the \code{plot3d.prcurve()} method uses the underlying
      functionality of the \pkg{rgl} package to provide an interactive
      3D representation of the data and the fitted principal curve.

      \code{prcurve()} also returns the fitted ordination and the
      original data, making it much easier to predict where passive
      samples should lay and simpler to plot the curve.

      \item New function (and S3 generic), \code{sppResponse()} for
      species responses along gradients. Currently the only supplied
      method is for the results of \code{prcurve()} fits, where the
      function returns the fitted response curves along the principal
      curve. A \code{plot()} method is also available.

      \item \code{logitreg()} can now fit the logistic regression via
      Firth's bias reduction method. This helps in cases where there is
      complete or quasi-complete separation, and for cases where there
      are far more no-analogues than analogues.

      \item This release sees much improved functionality for fitting
      principal component regression calibration models where the data
      are subject to an "ecologically meaningful" transformation (such
      as the Hellinger transformation). Models fitted with \code{pcr()}
      can be cross-validated using the new \code{crossval()} method, and
      predictions from the fitted model can be computed using the
      \code{predict()} method. Both ensure that the transformation is
      applied to the test or new data samples in the appropriate manner.

      \item New function \code{rankDC()} to compute the rank correlation
      between gradient distances (e.g. environmental variables) and
      distances in species composition. Has both base and Lattice
      graphics plot methods (the latter via \code{dotplot()}).

      \item \code{Stratiplot()} gains two new arguments: \code{labelAt}
      and \code{labelRot}, which allows control of the placement and
      rotation of variable labels on the panels of the plot when not
      using the strip.

      \item The \code{plot()} method for \code{timetrack()} now allows
      the plotting of linear combination or weighted average site
      scores.

      \item The \code{fitted()} method for \code{timetrack()} now has
      argument \code{which}, replacing the original \code{type}.

      \item There is a new \code{scores()} method for objects eturned by
      \code{timetrack()}.

      \item \code{analog()} gains a method allowing it to work directly
      on objects of class \code{"distance"}, as returned by
      \code{distance()}. This avoids recomputing the dissimilarity
      matrix if a suitable one is already available.
    }

  } % new features

  \subsection{Bug Fixes}{

    \itemize{
      \item \code{chooseTaxa()} no longer drops empty dimensions if the
      conditions supplied result in a single species being
      returned. Reported by Michael Burstert.

      A warning about presence of \code{NA}s when \code{na.rm = TRUE}
      was supplied is also fixed.

      \item The \code{formula} method of \code{timetrack()} was not
      well-implemented for various reasons, and could result in the
      species data not being transformed and other mistakes.

      The new code takes a one-sided formula describing the constraints
      for the underlying ordination, and variables in this formula will
      be searched for in the data frame supplied as argument \code{env}.

      \item \code{predict.wa()} would sometimes produce a 1-column
      matrix during the deshrinking step, which would result in an
      error. This empty dimension is now dropped so that a vector is
      always returned from the deshrinking function.
  
    }
    
  } % deprecated

  %% \subsection{DEPRECATED}{

  %%   \itemize{
  %%     \item
  
  %%   }
    
  %% } % deprecated
}