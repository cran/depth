\name{scontour}
\alias{scontour}
\alias{STD}
\alias{contourc}
\concept{Tukey spherical depth}
\title{Plotting spherical depth contours}
\description{
Traces spherical depth contours of a multivariate data set. Supports data on the circle or on the sphere.
}
\usage{
scontour(P, tracepoints=FALSE, colpoints="black", tracemed=TRUE,
   maxdepth=FALSE, xlim=c(0,2*pi), displaymed=FALSE,
   title="Circular Tukey contours", ylab="Tukey's circular depth",
   xlab=expression(theta), colmed=2, colarc="red", sizepoints=3)
}
\arguments{
  \item{P}{The data as a vector, a matrix, a data frame or a list.}
  \item{tracepoints}{Logical; if TRUE, data points are added to the plot.}
  \item{colpoints}{A specification for the color of the data points.}
  \item{tracemed}{Logical; if TRUE, the Tukey median is added to the plot.}
  \item{maxdepth}{Logical; On the circle only; if TRUE, the maximum depth is printed on the plot.}
  \item{xlim}{Numeric vectors of length 2, giving the \code{x} coordinate range.}
  \item{displaymed}{Logical; On the circle only; if TRUE, the median value is printed on the plot.}
  \item{title}{On the circle only, a specification for the plot title.}
  \item{ylab}{On the circle only, a specification for the y axis title.}
  \item{xlab}{On the circle only, a specification for the x axis title.}
  \item{colmed}{Color of the Tukey median on the plot.}
  \item{colarc}{On the sphere only, color of the spherical depth contours on the plot.}
  \item{sizepoints}{Size of plotted points.}
}
\details{Supports data on the circle or the sphere. For data on the circle, data must be expressed in polar coordinates as a angle in radians with values between 0 and \eqn{2\pi}. Data on the sphere can be expressed in Euclidean coordinates (\eqn{n} by 3 matrix) or in spherical coordinates (\eqn{n} by 2 matrix) where the first column contains \eqn{\theta} and the second column \eqn{\phi}.  The type of coordinates is determined automatically based on the dimensions of the input.}
\value{
\item{plot}{A plot of Tukey spherical depth if the input data are on the circle, or the Tukey spherical depth contours if the input data are on the sphere.}
If data are on the sphere only, a list of 3 elements is also outputted.
\item{1}{A sorted vector giving the depths of the plotted contours.}
\item{2}{A list of matrices with the vertices of every contour.}
\item{3}{The Euclidean coordinates of the Tukey median}
}
\references{
Liu, R.Y., Parelius, J.M. and Singh, K. (1999), Multivariate analysis by data depth: Descriptive statistics, graphics and inference (with discussion), \emph{Ann. Statist.}, \bold{27}, 783--858.

Mardia, K.V. and Jupp, E.J. (1999). \emph{Directional Statistics}, Wiley. 
}

\author{Maxime Genest.}

\seealso{\code{\link{sdepth}} for calculation of the depth of a point, \code{\link{smed}} for Tukey's spherical median.}
\examples{## Plot of Tukey spherical depth for data on the circle.
set.seed(2011)
scontour(runif(30,min=0,max=2*pi))

## Tukey spherical depth contours for data 
## on the shpere expressed in spherical coordinates.
scontour(cbind(runif(20,min=0,max=2*pi),runif(20,min=0,max=pi)))

## Tukey spherical depth contours for data 
## on the sphere expressed in Euclidean coordinates.
x=matrix(rnorm(60),ncol=3)
x=t(apply(x,1,function(y){y/sqrt(sum(y^2))}))
scontour(x)
}

\keyword{multivariate}
\keyword{nonparametric}
\keyword{robust}
\keyword{directional}
