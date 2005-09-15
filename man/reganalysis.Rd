\name{Regional Frequency Analysis}
\alias{reganalysis}
\title{Perform a whole Regional Frequency Analysis}
\description{Perform a whole Regional Frequency Analysis via the Index
  Flood model [Dalrymple, 1960]. The analysis is interactive and the
  end-user is asked for choices. It produces a test to detect if the
  region can been regarded as homogeneous or not, L-moment plots for
  graphical checks and estimation of the regional and at-site
  distribution.
}
\usage{reganalysis(data)}
\arguments{
  \item{data}{a list with each arguments correspond to of a site within
    the region plus argument 'record': a vector with the record for each
    site.}
}
\value{
  A whole regional frequency analysis is performed. That is a
  heterogeneity test for the selected region - see function
  \code{\link{heterogeneity}}, L-moments plots - see function
  \code{\link{lmomplots}}, and estimation of the regional and at-site
  distribution with eventually return level plots - see function
  \code{\link{regdist}, \link{locdist}}.
}
\examples{
data(region)
\dontrun{reganalysis(region)}
}
\author{Mathieu Ribatet}
\references{Dalrymple, T. (1960) \emph{Flood Frequency Analysis}. Water
  Supply.}
\keyword{htest}
    