\name{L-moments}
\alias{lmomreg}
\title{L-moments and statistics within a region.}
\description{Compute sample L-moments of each site within a region once
  normalized with a specified index flood, regional sample L-moments and
  3 statistics of interest.}
\usage{lmomreg(sample.sites,index.flood=mean)}
\arguments{
  \item{sample.sites}{A list giving data from each site.}
  \item{index.flood}{A function who computes the index flood from
    data. The default is to compute the sample mean.}
}
\value{
  A list with 3 arguments. \code{'site'}: a matrix containing the first
  4 sample L-moments of each site once normalized by the index
  flood, \code{'reg'}: the first 4 regional sample L-moments and
  \code{'V'}: 3 statistics of interest. For more information, take a
  look at the references.
}
\references{Hosking, J. R. M. and Wallis, J. R. (1997) \emph{Regional
    Frequency Analysis}. Cambridge University Press.
}
\author{Ribatet Mathieu}
\examples{
data(maxWind)
\dontrun{lmomreg(sample.site = maxWind)}
}
\keyword{htest}
