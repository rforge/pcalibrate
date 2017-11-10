\name{pCalibrate-package}
\alias{pCalibrate-package}
%\alias{pCalibrate}
\docType{package}
\title{
\packageTitle{pCalibrate}
}
\description{Implements transformations of one- and two-sided P-values
to minimum Bayes factors.
The minimum Bayes factor is the smallest possible Bayes factor
for the point null hypothesis against the alternative
within the specified class of alternatives.

The function \code{pCalibrate()} provides transformations (called calibrations in the sequel) 
for two-sided P-values
which consider the P-value as the data and are directly based on the distribution of the P-value 
under the null hypothesis and the alternative.
For one- and two-sided P-values from z-tests, \code{zCalibrate()} implements
calibrations for different classes of alternatives.
The function \code{tCalibrate()} provides the same functionality
for one- and two-sided P-values from t-tests.
The functions \code{FCalibrate()} and \code{LRCalibrate()} calibrate two-sided
P-values from the F-test of overall significance in the linear model or
likelihood ratio tests, respectively.
%\packageDescription{pCalibrate}
}
\details{
Package: pCalibrate

Type: Package

Title: Bayesian Calibrations of P-Values

Version: 0.1-2

Date: 2017-11-05

Author: Manuela Ott [aut, cre], Leonhard Held [aut] 

Maintainer: Manuela Ott <manuela.ott@uzh.ch>

Depends: epitools, exact2x2, MCMCpack

License: GPL (>=2)
%The DESCRIPTION file:
%\packageDESCRIPTION{pCalibrate}
%\packageIndices{pCalibrate}
%~~ An overview of how to use the package, including the most important ~~
%~~ functions ~~
}
\author{
Manuela Ott, Leonhard Held
%\packageAuthor{pCalibrate}
Maintainer: \packageMaintainer{pCalibrate}
}
\references{
% Berger J. O., Sellke T. (1987). Testing a point null hypothesis: The irreconcilability of P
% values and evidence (with discussion). \emph{Journal of the American Statistical Association},  \bold{82}, 112--139.
Held, L. and Ott, M. (2017). On p-values and Bayes factors.
 To appear in \emph{Annual Review of Statistics and Its Application, Volume 5}.
% ~~ Literature or other references for background information ~~
}
\keyword{ package
}
\keyword{htest}
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
pCalibrate(p=c(0.05, 0.01, 0.001), alternative="noninformative")
zCalibrate(p=c(0.05, 0.01, 0.005), type="one.sided", alternative="simple")
zCalibrate(p=c(0.05, 0.01, 0.005),  type="two.sided", alternative="normal")
tCalibrate(p=c(0.05, 0.01, 0.005), n=c(10, 20, 50), type="two.sided", alternative="normal")
FCalibrate(p=c(0.05, 0.01, 0.005), n=20, d=c(2, 5, 10), alternative="chi.squared")
LRCalibrate(p=c(0.05, 0.01, 0.005), df=2, alternative="simple")
}
