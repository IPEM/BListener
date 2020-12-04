\name{JustHockIt}
\alias{JustHockIt}
\docType{data}
\title{JustHockIt}
\description{
JustHockIt is the name of a dataset of duet singers who sang a composition in semi-hocket. Each duet performed eight times the same song. The dataset contains inter-onset-intervals and assessments of performance quality and agency.
}
\usage{data("JustHockIt")}
\format{
  A list with the following entries:
  \itemize{
    \item{\code{IOIdataset}, a list with numeric arrays of inter-onset-intervals in milliseconds for each performance}
    \item{\code{Condition}, a dataframe specifying whether the performance involved movement or no movement (= 2 conditions)}
    \item{\code{AllAnnotations}, a list with self-annotations of the performance quality, done by the performers themselves using a slider over time. These numbers are averages of the slider values}
    \item{\code{AllAgency}, a list containing numeric values representing the feeling of control (=agency) using a Likert scale}
}}
\details{
See \href{https://www.frontiersin.org/article/10.3389/fnins.2020.00619}{Dell' Anna et al. (2020)}.}
\references{
Leman, M.(2020).
Co-regulated timing in music ensembles - a Bayesian listener perspective. \emph{Submitted}.

Dell’Anna, A., Buhmann, J., Six, J., Maes, PJ., Leman, M., (2020).
\href{https://www.frontiersin.org/article/10.3389/fnins.2020.00619}{Timing Markers of Interaction Quality During Semi-Hocket Singing.} Frontiers in Neuroscience 14, 619
}

\keyword{datasets}
