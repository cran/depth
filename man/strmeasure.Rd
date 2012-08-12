\name{strmeasure}
\alias{strmeasure}
\concept{Tukey spherical depth}
\title{Computing trimmed measures of sherical location}
\description{Computes a sample trimmed measure of location based on the spherical Tukey's depth. Supports data on the circle or on the sphere (for Circular median only).}
\usage{
strmeasure(P,sorted=FALSE,depths=NULL,alpha=0,method="Mean")
}
\arguments{
  \item{P}{The data as a vector, a matrix, a data frame or a list. }
  \item{sorted}{Logical; if TRUE, it indicates that the data given in first argument is sorted.}
  \item{depths}{An optionnal vector of the same length of P that contains the Tukey's depth of each data. The calculation of the depth is then skipped and the provided values are used instead.}
  \item{alpha}{An optionnal numeric value between 0 and 1 to compute the median on a trimmed region rather than on the whole dataset. The trimming keeps only those points with a depth greater than or equal to \code{alpha}. The default value of \eqn{0} computes the median from the maximum depth trimmed region (i.e. no trimming).}
  \item{method}{Character string which determines the measure used. \code{method} can be "Mean" (the default) to compute trimmed mean direction or "Tukey" (for circular sample only) to compute trimmed Tukey's median.}
}
\details{This function returns a location estimate (Tukey's median or mean direction) of a sample truncated by Tukey's depth. For data on the circle, data must be expressed in polar coordinates as a angle in radians with values between 0 and \eqn{2\pi}. Data on the sphere can be expressed in Euclidean coordinates (\eqn{n} by 3 matrix) or in spherical coordinates (\eqn{n} by 2 matrix) where the first column contains \eqn{\theta} and the second column \eqn{\phi}. The type of coordinates is determined automatically based on the dimensions of the input.

While the option \code{method="Tukey"} supports only data on the circle, \code{method="Mean"} can also handle data on the sphere.}
\value{
   If the input sample is on the circle, a numeric value between \eqn{0} and \eqn{2\pi} giving the trimmed measure. If the input sample is on the sphere, the trimmed measure in Euclidean coordinates.
}

\references{Liu, R.Y., Parelius, J.M. and Singh, K. (1999), Multivariate analysis by data depth: Descriptive statistics, graphics and inference (with discussion), \emph{Ann. Statist.}, \bold{27}, 783--858.

Mardia, K.V. and Jupp, E.J. (1999). \emph{Directional Statistics}, Wiley. 
}

\author{Maxime Genest.}

\seealso{\code{\link{sdepth}} for the calculation of the depth of a point, \code{\link{scontour}} for Tukey's spherical median.}
\examples{## calculation of trimmed mean direction
set.seed(2011)
strmeasure(runif(30,min=0,max=2*pi),alpha=1/3,method="Mean")

## calculating of trimmed Tukey median
set.seed(2011)
strmeasure(runif(30,min=0,max=2*pi),alpha=1/3,method="Tukey")}

\keyword{multivariate}
\keyword{nonparametric}
\keyword{robust}
\keyword{directional}
