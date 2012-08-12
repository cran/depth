\name{sdepth}
\alias{sdepth}
\concept{Tukey spherical depth}
\title{Calculation of spherical depth}
\description{
Computes the spherical depth of a point with respect to a multivariate data set. Supports data on the circle or on the sphere.
}
\usage{
sdepth(theta, P)
}
\arguments{
  \item{theta}{Numerical vector whose depth is to be calculated. The coordinate system must match that of the observations.}
  \item{P}{The data as a vector, a matrix, a data frame or a list. }
}
\details{Computes the Tukey depth of \code{theta} with respect to the dataset \code{P}. For data on the circle, data must be expressed in polar coordinates as a angle in radians with values between 0 and \eqn{2\pi}. Data on the sphere can be expressed in Euclidean coordinates (\eqn{n} by 3 matrix) or in spherical coordinates (\eqn{n} by 2 matrix) where the first column contains \eqn{\theta} and the second column \eqn{\phi}. The type of coordinates is determined automatically based on the dimensions of the input.}
\value{
Returns the spherical depth of multivariate point \code{theta} with respect to the data set \code{P}.
}
\references{Liu, R.Y., Parelius, J.M. and Singh, K. (1999), Multivariate analysis by data depth: Descriptive statistics, graphics and inference (with discussion), \emph{Ann. Statist.}, \bold{27}, 783--858.

Mardia, K.V. and Jupp, E.J. (1999). \emph{Directional Statistics}, Wiley. 
}

\author{Maxime Genest.}

\seealso{\code{\link{scontour}} for depth graphics, \code{\link{smed}} for Tukey's spherical median.}
\examples{## Tukey spherical depth for a dataset on the circle
set.seed(2011)
sdepth(pi,runif(50,min=0,max=2*pi))

## Tukey spherical depth for data in spherical coordinates.
sdepth(c(pi,pi/2),cbind(runif(50,min=0,max=2*pi),runif(50,min=0,max=pi)))

## Tukey spherical depth for data in Eudlidean coordinates.
x=matrix(rnorm(150),ncol=3)
x=t(apply(x,1,function(y){y/sqrt(sum(y^2))}))
sdepth(x[1,],x)
}

\keyword{multivariate}
\keyword{nonparametric}
\keyword{robust}
\keyword{directional}
