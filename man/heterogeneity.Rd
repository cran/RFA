\name{Heterogeneity Statistics}
\alias{heterogeneity}
\title{Compute 3 heterogeneity statistics.}
\description{Compute 3 heterogeneity statistics based on the sample
  L-moments. It is usefull to test if a region could be considered as
  homogeneous or not, that is if all events from each site could came
  from only one distribution once normalized with the index flood.}
\usage{heterogeneity(N.sim,N.site,size,param,Vsite)}
\arguments{
  \item{N.sim}{Numeric. The number of stochastic region generated. It is
    recommended to set \code{N.sim} at list 500.}
  \item{N.site}{Numeric. The number of site within the region.}
  \item{size}{Numeric vector of length equal to \code{N.site}. The
    sample size of each site within the region.}
  \item{param}{Numeric vector of length 4. The parameter vector of the
    fitted regional Kappa distribution.}
  \item{Vsite}{Numeric vector of length 3. The statistic values for the
    region.}
}
\value{
  The program print value of the three statistics \code{H1}, \code{H2} and
  \code{H3} and return a list of several not so usefull objects.
}
\details{
  The code first fits a 4-parameter Kappa distribution from the weigthed
  mean of  the first 4 sample L-moments of each site. Then, it generates
  a large number of stochastics region with the same properties as the
  real one. Observations are sampled in the Kappa distribution fitted
  earlier.

  For each stochastic region, 3 statistics are computed and the
  behaviour of these values are directly compared with the 'real'
  statistics. Of course, if they have the same behaviour, then the
  region should be considered as homogeneous.

  The general form of the test statistics is \deqn{H = \frac{V -
      m}{\sigma}}{H = (V - m) / sigma} where \code{V} is one of the 3
  statistics \code{V1, V2, V3} defined later evaluated of the 'real'
  region, \code{m} and \eqn{\sigma}{sigma} are the mean and standart
  deviation of the \code{N} statistics \code{V} evaluated on each
  stochastic region.

  The 3 statistics \code{V1, V2} and \code{V3} are:
  \deqn{V1 = \left( \frac{\sum_{i=1}^N n_i \left(\tau_i - \tau_R
      \right)^2}{\sum_{i=1}^N n_i} \right)^{1/2}}{V = { \sum_{i=1}^N n_i
      ( tau_i - tau_R )^2 / \sum_{i=1}^N n_i }^(1/2)}
  \deqn{V2 = \frac{\sum_{i=1}^N n_i \left\{ \left(\tau_i - \tau_R \right)^2 +
    \left(\tau_{3,i} - \tau_{3,R} \right)^2 \right\}^{1/2}}{\sum_{i=1}^N
    n_i}}{\sum_{i=1}^N n_i { (tau_i - tau_R)^2 + (tau_(3,i) -
    tau_(3,R))^2 }^(1/2) / \sum_{i=1}^N n_i}
\deqn{V3 = \frac{\sum_{i=1}^N n_i \left\{ \left(\tau_{4,i} -
    \tau_{4,R} \right)^2 + \left(\tau_{3,i} - \tau_{3,R} \right)^2
    \right\}^{1/2}}{\sum_{i=1}^N n_i}}{\sum_{i=1}^N n_i { (tau_(4,i) -
    tau_(4,R))^2 + (tau_(3,i) - tau_(3,R))^2 }^(1/2) / \sum_{i=1}^N n_i}
  where \eqn{n_i}{n_i} is the sample size of site \code{i}, \eqn{\tau_i,
    tau_R}{tau_i, tau_R} are the at-site i and regional sample
  L-CV respectively,  \eqn{\tau_{3,i}, \tau_{3,R}}{tau_{3,i}, tau_{3,R}}
  are the at-site i and regional sample L-Skewness respectively,
  \eqn{\tau_{4,i}, \tau_{4,R}}{tau_{4,i}, tau_{4,R}} are the at-site
  i and regional sample L-Kurtosis respectively and N is
  the number of stochastic region generated.

  The last 2 statistics lack power to discriminate between homogeneous
  and heterogeneous region. They rarely yield values larger than 2 even
  for grossly heterogeneous regions.

  For the first statistic, it is suggest that the region be regarded as
  \emph{acceptably homogeneous} if \code{H1} < 1, \emph{possibly
    heterogeneous} if 1 \eqn{\le}{<=} \code{H1} < 2, \emph{definitely
    heterogeneous} if \code{H1} \eqn{\ge}{<=} 2.
}
\references{Hosking, J. R. M. and Wallis, J. R. (1997) \emph{Regional
    Frequency Analysis}. Cambridge University Press.
}
\author{Ribatet Mathieu}
\examples{
data(northCascades)
lmom <- c(1, 0.1103, 0.0279, 0.1366)
kappaParam <- kappalmom(lmom)
heterogeneity(500, 19, size = northCascades[,1],
kappaParam, c(0.0104, .0339, .0405))
##The heterogeneity statistics given by Hosking for this case
##study are H1 = 0.62, H2 = -1.49 and H3 = -2.37
##Taking into account sample variability, results should be
##consistent
}
\keyword{htest}
  