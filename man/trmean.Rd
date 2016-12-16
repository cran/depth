\name{trmean}
\alias{trmean}
\title{Classical-like depth-based trimmed mean}
\description{
Computes a sample trimmed mean based on the Tukey depth, the Liu depth or the Oja depth.
}
\usage{
trmean(x, alpha, W = function(dep, alpha){return(1)}, 
   method = "Tukey", ndir = 1000, approx = FALSE, 
   eps = 1e-8, ...)
}
\arguments{
  \item{x}{The data as a matrix, data frame or list. If it is a matrix or data frame, then each row is viewed as one bivariate observation. If it is a list, all components must be numerical vectors of equal length (coordinates of observations).}
  \item{alpha}{Outer trimming fraction (0 to 0.5). Observations whose depth is less than \code{alpha} to be trimmed.}
  \item{W}{Nonnegative weight function defined on [0, 1] through its argument \code{dep}. Number of arguments can be greater than 2 but the trimming fraction has to be one argument. See examples.}
  \item{method}{Character string which determines the depth function used. \code{method} can be "Tukey" (the default), "Liu" or "Oja".}
  \item{ndir}{Positive integer. Number of random directions used when approximate Tukey depth is utilised. Used jointly with \code{approx = TRUE}.}
  \item{approx}{Logical. If dimension is 3, should approximate Tukey depth be used? Useful when sample size is large.}
  \item{eps}{Error tolerance to control the calculation.}
  \item{...}{Any additional arguments to the weight function.}
}
\details{Dimension 2 or higher when \code{method} is "Tukey" or "Oja"; dimension 2 only when \code{method} is "Liu". Exactness of calculation depends on \code{method}. See \code{\link{depth}}. 
}
\value{
Multivariate depth-based trimmed mean}
\references{Masse, J.C and Plante, J.F. (2003), A Monte Carlo study of the accuracy and robustness of ten bivariate location estimators, \emph{Comput. Statist. Data Anal.}, \bold{42}, 1--26.

Masse, J.C. (2008), Multivariate Trimmed means based on the Tukey depth, \emph{J. Statist. Plann. Inference}, in press.

Rousseeuw, P.J. and Ruts, I. (1996), Algorithm AS 307: Bivariate location depth, \emph{Appl. Stat.-J. Roy. St. C}, \bold{45}, 516--526.}

\author{Jean-Claude Masse and Jean-Francois Plante, based on Fortran code by Ruts and Rousseeuw from University of Antwerp.}

\seealso{\code{\link{med}} for medians and \code{\link{ctrmean}} for a centroid trimmed mean.}
\examples{
## exact trimmed mean with default constant weight function
data(starsCYG, package = "robustbase")
trmean(starsCYG, .1)

## another example with default constant weight function
set.seed(159); library(MASS)
mu1 <- c(0,0); mu2 <- c(6,0); sigma <- matrix(c(1,0,0,1), nc = 2)
mixbivnorm <- rbind(mvrnorm(80, mu1, sigma), mvrnorm(20, mu2, sigma))
trmean(mixbivnorm, 0.3)

## trimmed mean with a non constant weight function
W1 <-function(x,alpha,epsilon) {
  (2*(x-alpha)^2/epsilon^2)*(alpha<=x)*(x<alpha+epsilon/2)+
  (-2*(x-alpha)^2/epsilon^2+4*(x-alpha)/epsilon-1)*
  (alpha+epsilon/2<=x)*(x<alpha+epsilon)+(alpha+epsilon<=x)
}
set.seed(345)
x <- matrix(rnorm(210), nc = 3)
trmean(x, .1, W = W1, epsilon = .05)

## two other examples of weighted trimmed mean
set.seed(345)
x <- matrix(rnorm(210), nc = 3)
W2 <- function(x, alpha) {x^(.25)}
trmean(x, .1, W = W2)
W3 <- function(x, alpha, beta){1-sqrt(x)+x^2/beta}
trmean(x, .1, W = W3, beta = 1)
}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{robust}
