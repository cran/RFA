\name{L-moments}
\alias{lmomplots}
\title{Generate several L-moments plots}
\description{Produces several plot given the first 3 sample L-moments
  ratio.}
\usage{lmomplots(lmom.reg, which = 1:3, ask = nb.fig < length(which) &&
dev.interactive(), pch = 15, draw.dist = TRUE)}
\arguments{
  \item{lmom.reg}{A list. Typically, it is the result of function
    \code{\link{lmomreg}}.}
  \item{which}{A numeric vector which specifies which L-moments plot
    should be plotted. '1' for L-Skewness - L-CV, '2' for L-Kurtosis -
    L-CV and '3' for L-Skewness - L-Kurtosis plot. By default, the 3
    plots are produced.}
  \item{ask}{Logical. If \code{TRUE} user is asked before a
    new plot is drawn.}
  \item{pch}{Numeric. Specifies the symbol to use for the 'regional
    point'. See \code{\link{points}}.}
  \item{draw.dist}{Logical. If set to \code{TRUE} - the default, then
    theoretical lines of some distribution are displayed in the
    L-Skewness / L-Kurtosis space.}
}
\value{
  Produces several plots.
}
\details{
  For the L-Skewness / L-Kurtosis plot, theoretical lines of Generalized
  Logistic, Generalized Extreme Value, Generalized Pareto, Log-Normal
  a,d Pearson Type III are displayed. OLB is the overall lower bound of
  \eqn{\tau_4}{tau_4} as a function of \eqn{\tau_3}{tau_3}. The
  triangles U, N, G and E correspond respectively to the Uniform,
  Normal, Gumbel and Exponential cases.
}
\author{Ribatet Mathieu}
\keyword{hplot}