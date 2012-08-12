\name{perspdepth}
\alias{perspdepth}
\title{ Perspective plots for depth functions }
\description{
Draws a perspective plot of the surface of a depth function over the x-y plane.
}
\usage{
perspdepth(x, method = "Tukey", output = FALSE, tt = 50, 
   xlab = "X", ylab = "Y", zlab = NULL, col = NULL, ...)
}
\arguments{
  \item{x}{Bivariate data as a matrix, data frame or list. If it is a matrix or data frame, then each row is viewed as one bivariate observation. If it is a list, both components must be numerical vectors of equal length (coordinates of observations).}
  \item{method}{Character string which determines the depth function used. \code{method} can be "Tukey" (the default), "Liu" or "Oja".}
  \item{output}{Logical. Default FALSE produces a perspective plot; otherwise, returns a list containing the grid points and depth values over these points.}
  \item{tt}{Gridsize. Number of equally spaced grid points in each coordinate direction to be used in perspective plot.}
  \item{xlab}{Title for x-axis. Must be a character string. }
  \item{ylab}{Title for y-axis. Must be a character string.}
  \item{zlab}{Title for z-axis. Must be a character string. Default NULL identifies the depth function.}
  \item{col}{ Color of the surface plot. Default NULL is "lightblue".}
  \item{...}{ Any additional graphical parameters.}
}
\details{Requires the \code{rgl} package. The perspective plot takes advantage of some of  the user interaction facilities of that package.
}
\value{Default \code{output = FALSE} yields a perspective plot; otherwise the function returns a list with components
  \item{x}{x-coordinates of the grid where the depth function is evaluated.}
  \item{y}{y-coordinates of the grid where the depth function is evaluated.}
  \item{z}{Matrix whose entry z[i ,j] is
 the value of the depth function at (x[i], y[j]).}
}
\references{Rousseuw, P.J. and Ruts, I. (1996), AS 307 : Bivariate location depth, \emph{Appl. Stat.-J. Roy. S. C}, \bold{45}, 516--526.}

\author{Jean-Claude Masse and Jean-Francois Plante, based on Fortran code by Rousseeuw, Ruts and Struyf from University of Antwerp.}

\seealso{\code{\link{isodepth}}, \code{\link{depth}} }
\examples{
## 2 perspective plots
data(geyser, package = "MASS")
perspdepth(geyser, col = "magenta")
set.seed(159); library(MASS)
mu1 <- c(0,0); mu2 <- c(6,0); sigma <- matrix(c(1,0,0,1), nc = 2)
mixbivnorm <- rbind(mvrnorm(80, mu1, sigma),mvrnorm(20, mu2, sigma))
perspdepth(mixbivnorm, col = "chartreuse")

## grid coordinates and corresponding depth values
set.seed(601)
x <- matrix(rnorm(48), nc = 2)
perspdepth(x, output = TRUE, tt = 10)
}

\keyword{multivariate}
\keyword{nonparametric}
\keyword{robust}
