\name{Simulation of random number}
\alias{rkappa}
\alias{rgpd}
\alias{rgev}
\alias{qgpd}
\alias{kappalmom}
\title{Generate random number and fit the 4 parameter Kappa distribution
  with L-moments.}
\description{Generate random number of the 4 parameter Kappa
  distribution. This distribution encompass the Generalized Extreme
  Value, Generalized Pareto and Generalized Logistic distributions.}
\usage{rkappa(n, loc, scale, shape1, shape2)
rgpd(n, loc = 0, scale = 1, shape = 0)
rgev(n, loc = 0, scale = 1, shape = 0)
qgpd(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)
kappalmom(lmom)
}
\arguments{
  \item{n}{Numeric. The number of pseudo-random number to be generated.}
  \item{p}{Numeric. The probability for which the quantile has to be
    computed.}
  \item{loc, scale, shape, shape1, shape2}{Numerics. The location, scale
    and shape parameters.}
  \item{lower.tail}{Logical. If \code{TRUE} the default, the quantile
    associated to probability \code{p} of non exceedance is
    computed. Otherwise, this is related to the probability of
    exceedance.}
  \item{lmom}{Numeric vector of length 4. The first 4 sample L-moments,
    that is the sample mean, L-CV, L-Skewness and L-Kurtosis.}
}
\value{
  The program return either a numeric vector containing estimates of the
  location, scale and shape parameter of the Kappa distribution either a
  fail flag informing a problem in the optimization algorithm.
  
  \code{'1'} L-moments invalid
  
  \code{'2'} \eqn{\left(\tau_3, \tau_4\right)}{(tau_3, tau_4)} lies
  above the generalized-logistic line (suggests that l-moments are not
  consistent with any kappa distribution with \code{shape2} > -1)
  
  \code{'3'} iteration failed to converge
  
  \code{'4'} unable to make progress from current point in iteration
  
  \code{'5'} iteration encountered numerical difficulties - overflow would
  have been likely to occur
  
  \code{'6'} iteration for the shape parameters converged, but overflow
  would have occurred when calculating location and scale parameters.
}
\details{
  The 4-parameter Kappa distribution has cumulative distribution function 
  \deqn{F\left(x\right) = \left[ 1 - shape2 \left(1 + \frac{shape1
      \left(x - loc \right)}{scale} \right)^{-1/shape1}
    \right]^{1/shape2}}{F(x) = [1 - shape2{1 + shape1(x - loc) / scale
    }^(-1/shape1)]^(1/shape2)}
  For \code{shape2} = -1, this is the Generalized Logistic distribution,
  for \code{shape2} = 0, the Generalized Extreme Value distribution and
  for \code{shape2} = 1, the Generalized Pareto distribution.
  
  Function \code{\link{kappalmom}} uses sample L-moments to fit the
  Kappa distribution. Newtow-Raphson iteration is used to solve the
  equations that express \eqn{\tau_3}{tau_3} and \eqn{\tau_4}{tau_4} as
  functions of \code{shape1} et \code{shape2}. \code{loc} and
  \code{scale} are calculated as functions of \eqn{\lambda_1,
    \tau,}{lambda_1, tau,} \code{shape1} and \code{shape2}. To ensure a
  1-1 relationship between parameters and L-moments, the parameter space
  is restricted. See reference for more details.
}
\references{Hosking, J. R. M. (1994) \emph{The four-parameter kappa
    distribution.} IBM Journal of Research and Development, \bold{38},
  251-258
}
\author{Ribatet Mathieu and Hosking J. R. M. for the fortran code.}
\examples{
lmom <- c(1, 0.1103, 0.0279, 0.1366)
kappalmom(lmom)
}
\keyword{distribution}