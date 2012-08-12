\name{depth-package}
\alias{depth-package}
\docType{package}
\title{
Depth functions tools for multivariate analysis
}
\description{This is a collection of functions applying depth functions methodology to multivariate analysis.  Besides allowing calculation of depth values and depth-based location estimators, the package includes functions for drawing contour plots and perspective plots of depth functions. 
}
\details{
\tabular{ll}{
Package: \tab depth\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2012-08-12\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
All functions apply to a multivariate data set. Function \code{\link{depth}} calculates the depth of a point with respect to the data set. Depth functions covered are Tukey's, Liu's and Oja's. Functions \code{\link{med}}, \code{\link{trmean}} and \code{\link{ctrmean}} return depth-based medians, classical-like trimmed means and centroid trimmed means, respectively. Functions \code{\link{perspdepth}} and \code{\link{isodepth}} draw perspective  and contour plots, respectively. Functions \code{\link{sdepth}}, \code{\link{smed}}, \code{\link{strmeasure}} and \code{\link{scontour}} give equivalent results for directional data.
}
\author{Maxime Genest <maxime.genest@clevislauzon.qc.ca>, Jean-Claude Masse <jcmasse@mat.ulaval.ca>, Jean-Francois Plante <jfplante@hec.ca>.

Maintainer: Jean-Francois Plante <jfplante@hec.ca>
}
\references{Liu, R.Y., Parelius, J.M. and Singh, K. (1999), Multivariate analysis by data depth: Descriptive statistics, graphics and inference (with discussion), \emph{Ann. Statist.}, \bold{27}, 783--858.

Liu, R.Y. and Singh, K. (1992), Directional data: Concepts of data depth on circles and spheres, \emph{Ann. Statist.}, \bold{20}, 1468--1484.

Mardia, K.V. and Jupp, E.J. (1999). \emph{Directional Statistics}, Wiley. 

Small, C.G. (1990), A survey of multidimensional medians, \emph{Int. Statist. Rev.}, \bold{58}, 263--277.

Zuo, Y. amd Serfling, R. (2000), General Notions of Statistical Depth Functions, \emph{Ann. Statist.}, \bold{28}, no. 2, 461--482.}

\keyword{package}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{robust}

%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
set.seed(159); library(MASS)
mu1 <- c(0,0); mu2 <- c(6,0); sigma <- matrix(c(1,0,0,1), nc = 2)
mixbivnorm <- rbind(mvrnorm(80, mu1, sigma), mvrnorm(20, mu2, sigma))
depth(c(0,0),mixbivnorm)
med(mixbivnorm)
trmean(mixbivnorm, 0.2)
library(rgl)
perspdepth(mixbivnorm, col = "magenta")
isodepth(mixbivnorm, dpth = c(35,5), col = rainbow(2))
}
