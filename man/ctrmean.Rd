\name{ctrmean}
\alias{ctrmean}
\alias{ctrmean.default}
\alias{ctrmean.data.frame}
\alias{ctrmean.list}
\title{Centroid trimmed mean}
\description{
Computes the centroid of a Tukey depth-based trimmed region.
}
\usage{
ctrmean(x ,alpha, eps = 1e-8, mustdith = FALSE, maxdith = 50, 
   dithfactor = 10 ,factor = .8)
}
\arguments{
  \item{x}{Bivariate data as a matrix, data frame or list. If it is a matrix or data frame, then each row is viewed as one bivariate observation. If it is a list, both components must be numerical vectors of equal length (coordinates of observations).}
 \item{alpha}{Outer trimming fraction (0 to 0.5). Observations whose depth is less than \code{alpha} to be trimmed.}
  \item{eps}{Error tolerance to control the calculation.}
  \item{mustdith}{Logical. Should dithering be applied? Used when data set is not in general position or a numerical problem is encountered.}
  \item{maxdith}{Positive integer. Maximum number of dithering steps.}
  \item{dithfactor}{Scaling factor used for horizontal and vertical dithering.}
\item{factor}{Proportion (0 to 1) of outermost contours computed according to a version of the algorithm ISODEPTH of Rousseeuw and Ruts (1998); remaining contours are derived from an algorithm in Rousseeuw \emph{et al.} (1999). }
}
\details{Dimension 2 only. Centroid trimmed mean is defined to be the centroid
of a Tukey depth-based  trimmed region relative to the uniform measure. Contours
are derived from algorithm ISODEPTH by Ruts and Rousseeuw (1996) or, more
exactly, revised versions of this algorithm which appear in Rousseeuw and Ruts
(1998) and Rousseeuw \emph{et al.} (1999). Argument \code{factor} determines
which version to use. If \eqn{n} is the number of observations, contours of
depth \eqn{\le } \code{factor} \eqn{n/2} are obtained from the 1998 version, while
the remaining contours are derived from the 1999 version. 

When the data set is not in general position, dithering can be used in the sense that random noise is added to each component of each observation. Random noise takes the form \code{eps} times \code{dithfactor} times U for the horizontal component and \code{eps} times \code{dithfactor} times V for the vertical component, where U, V are independent uniform on [-.5, 5.]. This is done in a number of consecutive steps applying independent U's and V's.}
\value{Centroid trimmed mean vector
}
\references{Masse, J.C. (2008), Multivariate Trimmed means based on the Tukey depth, \emph{J. Statist. Plann. Inference}, in press.

Ruts, I. and Rousseeuw, P.J. (1996), Computing depth contours of bivariate point clouds, \emph{Comput. Statist. Data Anal.}, \bold{23}. 153--168.

Rousseeuw, P.J. and Ruts, I. (1998), Constructing the bivariate
Tukey median, \emph{Stat. Sinica}, \bold{8}, 828--839.

Rousseeuw, P.J., Ruts, I., and Tukey, J.W. (1999), The Bagplot: A Bivariate Boxplot, \emph{The Am. Stat.}, \bold{53}, 382--387.}

\author{Jean-Claude Masse and Jean-Francois Plante, based on Fortran code by Ruts and Rousseeuw from University of Antwerp.}

\seealso{\code{\link{med}} for multivariate medians and \code{\link{trmean}} for classical-like depth-based trimmed means.}
\examples{
## exact centroid trimmed mean
set.seed(345)
xx <- matrix(rnorm(1000), nc = 2)
ctrmean(xx, .2)

## second example of an exact centroid trimmed mean
set.seed(159); library(MASS)
mu1 <- c(0,0); mu2 <- c(6,0); sigma <- matrix(c(1,0,0,1), nc = 2)
mixbivnorm <- rbind(mvrnorm(80, mu1 ,sigma), mvrnorm(20, mu2, sigma))
ctrmean(mixbivnorm, 0.3)

## dithering used for data set not in general position
data(starsCYG, package = "robustbase")
ctrmean(starsCYG, .1, mustdith = TRUE)
}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{robust}
