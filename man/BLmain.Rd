\name{BLmain}
\alias{BLmain}
\title{Bayesian listener main function}
\description{The \code{BLmain} function does multivariate tracking of timing in inter-onset intervals.}
\usage{
BLmain(u, meter, tg = NULL, outt = 0.7, ssa = 100,
       g = 0, V = 0.00001, W = 0.01)
}
\arguments{
    \item{u}{a numeric vector indicating IOIs in milliseconds}
    \item{meter}{a numeric vector indicating the expected out-of-time IOI-meter}
    \item{tg}{a number indicating the expected duration of the smallest IOl-class in milliseconds}
    \item{outt}{a number specifying the threshold for outliers (default: 1.5 log2dur)}
    \item{ssa}{a number indicating the samplingrate of the system (default: 100 sa/sec)}
    \item{g}{a number of the regularisation strength applied to IOI-classes (default: 0, max=1)}
    \item{V}{a number specifiing the uncertainty in the system dynamics (default: 0.00001)}
    \item{W}{a number specifying the uncertainty in the observation (default: 0.01)}
}
\details{
\code{BLmain} does multivariate tracking of timing, given fluctuation in inter-onset-intervals. The tracking is based on Bayesian inferencing on latent constancy values, which predict future, to be observed, inter-onset-intervals over time. Relevant concepts are:
\itemize{
\item \emph{IOI-observation}, or \emph{IOI-observed object}: The time between two successive onsets.
\item \emph{IOI-class}: Constancy extracted from the IOI-observations using Bayesian inference. \emph{Multiple IOI-classes} can be tracked in parallel, resulting in a multivariate time series of timing constancy.
\item \emph{Regularization}: Bayesian inferencing on  IOI-classes can be regularized using a hyper-parameter to constraints variance among IOI-classes. Regularization uses the \emph{out-of-time IOI-meter}, which defines an assumed ratio among IOI-classes, disregarding how time evolves.
\item \emph{Time-related IOI-meter}: The ratios of the multivariate IOI-class over time.
\item \emph{IOI-tempo}: A linear combination of the IOI-classes at any moment in time, using the out-of-time IOI-meter ratios to calculate a mean of the IOI-class time series.
}
The algorithm obeys dynamic laws with updating at a samplingrate of 100 samples per second, unless otherwise stipulated (up to 1000 samples per second is possible).  However, IOI objects are handled in log2dur units, a log2 transformation of the millisecond scale, with "dur" standing for "duration".

}
\value{
a list with the following named entries:
    \item{Observed}{a matrix with IOI-observed objects}
    \item{Innovation}{a matrix with prediction errors, i.e. difference of predicted IOI based on the IOI-class, and the IOI-observed}
    \item{AdaptiveGain}{a matrix specifying the strength of adaptation (also: Kalman gain)}
    \item{Prediction}{a matrix with IOI-classes}
    \item{Tempo}{a numerical vector estimating tempo}
    \item{Outliers}{a matrix with outliers, or IOI-observations that fall outside the outlier threshold}
    \item{parameters}{a list of parameters given in the input function}
}
% It is better to define meter in ratio of milliseconds and do the log2 transformation inside the function.
\references{
Leman, M.(2020). Co-regulated timing in music ensembles - a Bayesian listener perspective.


\emph{Submitted}.
}
\author{Marc Leman}
\examples{
a <- c(205, 407, 395, 801,200, 412, 399, 810, 195,
       415, 400, 821,190, 420, 395, 835)
b <- c(4,2,1)
BLmain(a, b)}

\keyword{ Bayesian inference }
