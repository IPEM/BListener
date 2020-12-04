\name{BListener-package}
\alias{BListener-package}
\alias{BListener}
\docType{package}
\title{Analysis of co-regulated timing in music ensembles}
\description{Analysis of co-regulated timing in music ensembles}
\details{
\tabular{ll}{
Package: \tab BListener\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2020-11-20\cr
License: \tab GPL-2\cr
URL:
     \tab \url{https://github.com/ArtScienceLab/BListener} \cr
}

BListener is an algorithm that simulates how a listener-musician in a music ensemble captures the timing of that ensemble, such that actions could be taken to co-regulate the timing. BListener implements the perception part, not the action part. The algorithm is conceived to work in real-time although further work is needed to use it in a real-time application.

The input is an array of inter-onset intervals (called IOI-observations) whose cummulative sum corresponds with onset times in milliseconds. Those onsets can be obtained by onset-extractor algorithms, or by mannual annotation.

The main outputs are multivariate time series that predict the timing at each sample (IOI-classes), and the associated prediction error time series, which contain the difference between the predicted timing (IOI-class for time t) and the observed timing (IOI-observed objet at time t). From those time series it is possible to extract gobal features such as stability (standard deviation of IOI-classes over time) and fluctuation (mean of absolute values of prediction errors).

BListener is described in the article “Co-regulated timing in music ensembles - a Bayesian listener perspective”. To understand the terminology about IOI-objects, you should read this text. All functions and datasets are available. Figures and tables can be reproduced and we invite the reader to work out different applications.\cr
%example of cited function \code{\link{BLmain}}.
}
\author{Marc Leman; credits are also due to Marc Vidal.}
\references{
%Exaample with link
%Leman, M., Lesaffre, M., Tanghe, M. (2001).
%\href{https://biblio.ugent.be/publication/135924}{Introduction to the IPEM Toolbox for Perception-based Music Analysis.}
%\emph{Proceedings of the XIII Meeting of the FWO Research Society on Foundations of Music Research}, Ghent.

Leman, M.(2020).
Co-regulated timing in music ensembles - a Bayesian listener perspective. \emph{Submitted}.
}

\keyword{ Bayesian inference }
% \seealso{ }

