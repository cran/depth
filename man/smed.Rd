\name{smed}
\alias{smed}
\concept{Tukey spherical median}
\title{Calculating spherical medians}
\description{Computes the spherical median of a data set on the circle.}
\usage{
smed(P, sort=FALSE, depths=NULL, alpha=NULL,
   method="Tukey", tracecontour=FALSE, tracepoints=FALSE)
}
\arguments{
  \item{P}{The data as a vector, a matrix, a data frame or a list.}
  \item{sort}{Logical; TRUE indicates that the data in \code{P} is already sorted.} 
  \item{depths}{For Tukey's method only; An optionnal vector of the same length as P that contains the Tukey depth of each data. The calculation of the depth is then skipped and the provided values are used instead.} 
  \item{alpha}{For Tukey's method only; \code{alpha} is an optionnal numeric value between 0 and 1 to compute the median on a trimmed region rather than on the whole dataset. The trimming keeps only those points with a depth greater than or equal to \code{alpha}. The default value of \code{NULL} computes the median from the maximum depth trimmed region (i.e. no trimming).}
  \item{method}{Character string which determines the depth function used. \code{method} can be "Tukey" (the default) or "Circular".}
  \item{tracecontour}{Only if \code{method="Circular"}. Traces the plot of depth with respect to angular positions on the circle.}
  \item{tracepoints}{Only if \code{method="Circular"}. Draws the points and their median on the circle.}
}
\details{Calculates shperical medians for data on the circle only. The input must be a list of angles in radians between 0 and \eqn{2\pi} (polar coordinates). If \code{method="Tukey"}, the Tukey median is returned. If \code{method="Circular"}, the circular median (the point minimizing the average distance based on arccosine) is returned.}
\value{
A numeric value between \eqn{0} and \eqn{2\pi} giving the median in polar coordinate.
}
\references{Liu, R.Y., Parelius, J.M. and Singh, K. (1999), Multivariate analysis by data depth: Descriptive statistics, graphics and inference (with discussion), \emph{Ann. Statist.}, \bold{27}, 783--858.

Mardia, K.V. and Jupp, E.J. (1999). \emph{Directional Statistics}, Wiley. 
}

\author{Maxime Genest.}

\seealso{\code{\link{sdepth}} for calculation of the depth of a point, \code{\link{scontour}} for Tukey's spherical median.}
\examples{## calculation of the Tukey spherical median for data on the circle
set.seed(2011)
smed(runif(30,min=0,max=2*pi))
}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{robust}
\keyword{directional}
