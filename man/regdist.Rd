\name{Regional Frequency Analysis}
\alias{regdist}
\alias{locdist}
\title{Regional and at-site distribution}
\description{The regional and at-site distribution are estimated via the
  Index Flood model. This is achieved with the regional L-moments.}
\usage{regdist(lmom, loc, main, xlab, ylab, xliminf, xlimsup, draw.site =
  TRUE, \dots)
  locdist(param.reg, mu, data, main, xlab, ylab, xlimsup, draw.data = TRUE,
  \dots)
}
\arguments{
  \item{lmom}{A list. Typically, return of function
    \code{\link{lmomreg}}.}
  \item{param.reg}{Vector of length 3 giving the value of the parameter
    of the regional distribution.}
  \item{mu}{Numeric. The mean number of events per year.}
  \item{loc}{Optional numeric. If present, the location parameter is
    fixed to this value.}
  \item{data}{Vector of the at-site sample.}
  \item{main, xlab, ylab, xliminf, xlimsup}{Optional. Specifies the title, labels
    for the x and y axis and the lower and upper bound for the x-axis.}
  \item{index.flood}{A function who compute the index flood in functions
    of the at-site sample. The default is to compute the sample mean.}
  \item{draw.site, draw.data}{Logical. Should adimensional distribution
    of all site and at-site sample be plotted. The default is to plot
    it.}
  \item{\dots}{Optional parameters to be passed to the
    \code{\link{plot}} function.}
}
\value{
  The program return the parameter estimate of the regional and at-site
  distribution. Optionally, return level plots are displayed.} 
\author{Mathieu Ribatet}
\details{
  The Index Flood model [Dalrymple, 1960] states that all distribution within a
  homogeneous region are the same up to a normalisation by a constant -
  namely, the Index Flood. So, estimate the distribution of a specified
  site can be performed in two step:
  \itemize{
    \item Estimation of the regional distribution.
    \item Estimation of the at-site distribution by multiplying the regional one with the at-site index flood.
  }

  The regional distribution is estimated with the mean of L-moments of
  each site weighted with the record length [Hosking and Wallis, 1997].
}
\references{
  Dalrymple, T. (1960) \emph{Flood Frequency Analysis}. Water Supply
  Paper.
  
  Hosking, J. R. M. and Wallis, J. R. (1997) \emph{Regional Frequency
    Analysis}. Cambridge University Press.
}
\keyword{htest}
