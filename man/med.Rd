\name{med}
\alias{med}
\alias{med.default}
\alias{med.data.frame}
\alias{med.list}
\title{Multivariate median}
\description{
Computes the median of a multivariate data set.
}
\usage{
med(x, method = "Tukey", approx = FALSE, eps = 1e-8, maxit = 200, mustdith = FALSE, maxdith = 50, dithfactor = 10, factor = 0.8, nstp = floor(5*n^0.3*p), ntry = 10*(p+1), nalt = 4*(p+1), ndir = 1000)
}
\arguments{
  \item{x}{The data as a matrix, data frame or list. If it is a matrix or data frame, then each row is viewed as one multivariate observation. If it is a list, all components must be numerical vectors of equal length (coordinates of observations).}
  \item{method}{Character string which determines the depth function used. \code{method} can be "Tukey" (the default), "Liu", "Oja", "Spatial" or "CWmed".}
  \item{approx}{Logical. Should an approximate Tukey median be computed? Useful in dimension 2 only when sample size is large.}
  \item{eps}{Error tolerance to control the calculation.}
  \item{maxit}{Number of Newton-Raphson iterations in case \code{method} is "Spatial".}
  \item{mustdith}{Logical.Should dithering be applied? Used to compute the Tukey median when data set is not in general position or a numerical problem is encountered.}
  \item{maxdith}{Integer. Maximum number of dithering steps.}
  \item{dithfactor}{Scaling factor used for horizontal and vertical dithering.}
 \item{factor}{Proportion (0 to 1) of outermost contours computed according to algorithm HALFMED of Rousseeuw and Ruts (1998); remaining contours derived from an algorithm in Rousseeuw \emph{et al.} (1999).}
  \item{nstp}{Positive integer. Maximum number of steps in the iteration process leading to an approximate value of the Tukey median. Default value is taken to be the largest integer not greater than \eqn{5 n^{0.3}p}, where \eqn{n} is the number of observations and \eqn{p} the dimension.}
  \item{ntry}{Positive integer. Maximum number of steps without an increase of the Tukey depth in the iteration process leading to an approximate value of the Tukey median. Default value is taken to be \eqn{10(p+1)}, where \eqn{p} is the dimension.}
  \item{nalt}{Positive integer. Maximum number of consecutive steps without an increase of the Tukey depth at any time in the iteration process leading to an approximate value of the Tukey median. Default value is taken to be \eqn{4(p+1)}, where \eqn{p} is the dimension.}
  \item{ndir}{Positive integer. Number of random directions used in the iteration process leading to an approximate value of the Tukey median.}
}
\details{\code{method} "Tukey" computes the Tukey median. Calculation is exact
in dimensions 1 and 2, and approximate in higher dimensions. The bivariate case
utilises algorithm HALFMED by Rousseeuw and Ruts (1998) as well as an algorithm
from Rousseeuw \emph{et al.} (1999). Argument \code{factor} determines which
algorithm to use. If \eqn{n} is the number of observations,  contours of depth
\eqn{\le} \code{factor} \eqn{n/2} are derived from algorithm HALFMED, while the remaining contours are obtained from the second algorithm. The higher dimensional case is covered by Fortran code from Struyf and Rousseeuw (2000).

When \code{method} is "Tukey", data must be in general position. If not, in dimension 2 dithering can be used in the sense that random noise is added to each component of each observation. Random noise takes the form \code{eps} times \code{dithfactor} times U for the horizontal component and \code{eps} times \code{dithfactor} times V for the vertical component, where U, V are independent uniform on [-.5, 5.]. This is done in a number of consecutive steps applying independent U's and V's.

\code{method} "Liu" computes the Liu median. It is based on Fortran code from Rousseeuw and Ruts (1996) and restricted to two-dimensional data.

\code{method} "Oja" computes the Oja median. It is based on Fortran code by Niinimaa et al. (1992) and restricted to two-dimensional data.

\code{method} "Spatial" computes the spatial median or mediancentre. It is based on Fortran code by Gower (1974), and Bedall and Zimmermann (1979). 

\code{method} "CWmed" computes the coordinatewise median.
}
\value{A list with components
  \item{median}{the median}
  \item{depth}{the depth of the median (omitted when \code{method} is "Spatial" or "CWmed")}
}
\references{Gower, J.C. (1974), AS 78: The Mediancentre, \emph{Appl. Stat.}, \bold{23}, 466--470. 

Bedall, F.K. and Zimmermann, H. (1979), AS 143: The Mediancentre, \emph{Appl. Stat.}, \bold{28}, 325--328.

Niinimaa, A, Oja, H., Nyblom, J. (1992), AS 277 : The Oja Bivariate Median, \emph{Appl. Stat.}, \bold{41}, 611--617.

Rousseeuw, P.J. and Ruts, I. (1996), Algorithm AS 307: Bivariate location depth, \emph{Appl. Stat.-J. Roy. St. C}, \bold{45}, 516--526.

Rousseeuw, P.J. and Ruts, I. (1998), Constructing the bivariate
Tukey median, \emph{Stat. Sinica}, \bold{8}, 828--839.

Rousseeuw, P.J., Ruts, I., and Tukey, J.W. (1999), The Bagplot: A Bivariate Boxplot, \emph{The Am. Stat.}, \bold{53}, 382--387.

Small, C.G. (1990), A survey of multidimensional medians, \emph{Int. Statist. Rev.}, \bold{58}, 263--277.

Struyf, A. and Rousseeuw, P.J. (2000), High-dimensional computation of the deepest location, \emph{Comput. Statist. Data Anal.}, \bold{34}, 415--436.

Masse, J.C and Plante, J.F. (2003), A Monte Carlo study of the accuracy and robustness of ten bivariate location estimators, \emph{Comput. Statist. Data Anal.}, \bold{42}, 1--26.
}
\author{Jean-Claude Masse and Jean-Francois Plante, based on Fortran code by authors listed in the references.}

\seealso{\code{\link{trmean}} and \code{\link{ctrmean}} for trimmed means}
\examples{
## exact Tukey median for a mixture of bivariate normals
set.seed(159); library(MASS)
mu1 <- c(0,0); mu2 <- c(6,0); sigma <- matrix(c(1,0,0,1), nc = 2)
mixbivnorm <- rbind(mvrnorm(80, mu1, sigma), mvrnorm(20, mu2, sigma))
med(mixbivnorm)

##  approximate Tukey median of a four-dimensional data set
set.seed(601)
zz <- matrix(rnorm(96), nc = 4)
med(zz)

## data set not in general position
data(starsCYG, package = "robustbase")
med(starsCYG, method = "Liu")

## use of dithering for the Tukey median
med(starsCYG, mustdith = TRUE)
}
\keyword{multivariate}
\keyword{nonparametric}
\keyword{robust}
