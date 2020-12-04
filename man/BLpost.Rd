\name{BLpost}
\alias{BLpost}
\title{Postprocessing of "BL" objects}
\description{The function \code{BLpost} extracts global measures of timing}
\usage{BLpost(R)}
\arguments{
  \item{R}{An object of class "BL"}
}
\details{
The function \code{BLpost} extracts global measures of timing, using a BL object as input; see \code{\link{BLmain}}.
The relevant concepts are fluctuation, stability, narration, collapse, and constancy.
\itemize{
\item \emph{Fluctuation} is for each IOI-class a global measure of its prediction error (cf. the Innovation matrix in the BL object). \emph{Fluctuation1} calculates the mean (over time) of the absolute Innovation values, whereas \emph{Fluctuation2} calculates their standard deviation (over time).
\item \emph{Stability} is for each IOI-class a global measure of the standard deviation of over time.
\item \emph{Collapse} is a global measure of outlier time, in milliseconds.
\item \emph{Constancy} is a global measure of variance of the time-related IOI-meter. This only applies to multivariate BL objects (minimum 2 IOI-classes). Under development}
}
\value{
a list with the following named entries:
    \item{Fluctuation1}{A vector representing means, using absolute values of the Innovation matrix}
    \item{Fluctuation2}{A vector representing standard deviations, using the Innovation matrix}
    \item{Stability}{A vector with standard deviations, using the Prediciton matrix}
    \item{Collapse}{A numerical value representing the total duration of outlier objects}
    \item{Constancy}{A numerical value representing the mutivariate variance in the Prediction matrix}
}
\references{
Leman, M.(2020). Co-regulated timing in music ensembles - a Bayesian listener perspective.
\emph{Submitted}.
}
\author{Marc Leman}
\examples{
a <- c(205, 407, 395, 801,200, 412, 399, 810, 195,
       415, 400, 821,190, 420, 395, 835)
b <- c(4,2,1)
BL.results <- BLmain(a, b)
post.R <- BLpost(BL.results)
}
