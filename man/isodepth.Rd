\name{isodepth}
\alias{isodepth}
\title{Contour plots for depth functions}
\description{
Draws a contour plot of Tukey's depth function.
}
\usage{
isodepth(x, dpth = NULL, output = FALSE, twodim = TRUE, 
   mustdith = FALSE, maxdith = 50, dithfactor = 10, 
   trace.errors = TRUE, eps = 1e-8, factor = 0.8, xlab = "X", 
   ylab = "Y", zlab = "Tukey's depth", colcontours = NULL, ...)
}
\arguments{
  \item{x}{Bivariate data as a matrix, data frame or list. If it is a matrix or data frame, then each row is viewed as one bivariate observation. If it is a list, both components must be numerical vectors of equal length (coordinates of observations).}
  \item{dpth}{Vector of positive integers. Numbers \eqn{1,2,\ldots} refer to contours of depth \eqn{1/n, 2/n,\ldots}, where \eqn{n} is the number of observations. Useful to draw particular contours. Default \code{dpth  = NULL} corresponds to the set of all contours.}
  \item{output}{Logical. Default FALSE produces a contour plot; otherwise a list of contour vertices. }
  \item{twodim}{Logical. \code{twodim = FALSE} returns a transparent perspective plot making use of the \code{rgl} package. }
  \item{mustdith}{Logical. Should dithering be applied? Used when data set is not in general position or a numerical problem is encountered.}
  \item{maxdith}{Positive integer. Maximum number of dithering steps.}
  \item{dithfactor}{Scaling factor used for horizontal and vertical dithering.}
  \item{trace.errors}{Logical. Should all contours be considered? Used when a numerical problem is encountered for some inner contours. Default \code{trace.errors = FALSE} means those contours are left out.}
  \item{eps}{Error tolerance to control the calculation. }
  \item{factor}{Proportion (0 to 1) of outermost contours computed according to a version of the algorithm ISODEPTH of Rousseeuw and Ruts (1998); remaining contours are derived from an algorithm in Rousseeuw \emph{et al.} (1999). }
  \item{xlab}{Title for x-axis. Must be a character string. }
  \item{ylab}{Title for y-axis. Must be a character string.}
  \item{zlab}{Title for z-axis. Used jointly with \code{twodim = FALSE}.}
  \item{colcontours}{Vector of color names of some or all of the contours. Recycling is used when necessary. Colors can be specified in different ways, see color specification in \code{\link{par}},}
  \item{...}{Any additional graphical parameters (see \code{par}). }
}
\details{
Tukey's depth and dimension 2 only. Contours are computed according to algorithm
ISODEPTH by Ruts and Rousseeuw (1996) or, more exactly, revised versions of this
algorithm which appear in Rousseeuw and Ruts (1998) and Rousseeuw \emph{et al.}
(1999). Argument \code{factor} determines which version to use. If \eqn{n} is
the number of observations, contours of depth \eqn{\le} \code{factor} \eqn{n/2} are obtained from the 1998 version, while the remaining contours are derived from the 1999 version.

When the data set is not in general position, dithering can be used in the sense that random noise is added to each component of each observation. Random noise takes the form \code{eps} times \code{dithfactor} times U for the horizontal component and \code{eps} times \code{dithfactor} times V for the vertical component, where U, V are independent uniform on [-.5, 5.]. This is done in a number of consecutive steps applying independent U's and V's.
}
\value{Default \code{output = FALSE} yields a contour plot. If not, the function returns a list of \eqn{m} components, where \eqn{m} is the number of contours and component \eqn{i} is a matrix whose rows are the vertices of contour \eqn{i}.
}
\references{Ruts, I. and Rousseeuw, P.J. (1996), Computing depth contours of bivariate point clouds, \emph{Comput. Stat. Data An.}, \bold{23}, 153--168.

Rousseeuw, P.J. and Ruts, I. (1998), Constructing the bivariate
Tukey median, \emph{Stat. Sinica}, \bold{8}, 828--839.

Rousseeuw, P.J., Ruts, I., and Tukey, J.W. (1999), The Bagplot: A Bivariate Boxplot, \emph{The Am. Stat.}, \bold{53}, 382--387.}

\author{Jean-Claude Masse and Jean-Francois Plante, based on Fortran code by Rousseeuw, Ruts and Struyf from University of Antwerp.}
\seealso{\code{\link{depth}}, \code{\link{perspdepth}} }
\examples{
## exact contour plot with 10 contours
set.seed(601) ; x = matrix(rnorm(48), nc = 2)
isodepth(x)

## exact colored contours
set.seed(159); library(MASS)
mu1 <- c(0,0); mu2 <- c(6,0); sigma <- matrix(c(1,0,0,1), nc = 2)
mixbivnorm <- rbind(mvrnorm(80, mu1 ,sigma), mvrnorm(20, mu2, sigma))
isodepth(mixbivnorm, dpth = c(35,5), col = rainbow(2))

## vertices of each contour
set.seed(601)
x <- matrix(rnorm(48), nc = 2)
isodepth(x, output = TRUE)

## data set not in general position
data(starsCYG, package = "robustbase")
isodepth(starsCYG, mustdith = TRUE)

## colored contours
set.seed(601)
x <- matrix(rnorm(48), nc = 2)
isodepth(x, colcontours= rainbow(10))

# perspective plot
library(rgl)
set.seed(601)
x <- matrix(rnorm(48), nc = 2)
isodepth(x, twodim = FALSE)
}

\keyword{multivariate}
\keyword{nonparametric}
\keyword{robust}
