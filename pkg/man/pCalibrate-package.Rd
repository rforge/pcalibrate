\name{pCalibrate-package}
\alias{pCalibrate-package}
%\alias{pCalibrate}
\docType{package}
\title{
\packageTitle{pCalibrate}
}
\description{Implements transformations of one- and two-sided p-values to minimum Bayes factors. The minimum Bayes factor is the smallest possible Bayes factor for the point null hypothesis against the alternative within the specified class of alternatives.

The function \code{pCalibrate()} provides minimum Bayes factors
for two-sided p-values
which consider the p-value as the data and are directly based on the distribution of the p-value 
under the null hypothesis and the alternative.
For one- and two-sided p-values from z-tests, \code{zCalibrate()} implements
minimum Bayes factors for different classes of alternatives.
The function \cr 
\code{tCalibrate()} provides the same functionality
for one- and two-sided p-values from t-tests.
The functions \code{FCalibrate()} 
and \code{LRCalibrate()} transform two-sided
p-values from the F-test or likelihood ratio test, respectively, to
minimum Bayes factors.
%\packageDescription{pCalibrate}
}
\details{
Package: pCalibrate

Type: Package

Title: Bayesian Calibrations of p-Values

Version: 0.2-1

Date: 2020-03-19

Author: Manuela Ott [aut, cre], Leonhard Held [aut] 

Maintainer: Manuela Ott <manuela.c.ott@gmail.com>

Depends: exact2x2, MCMCpack

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
Held, L. and Ott, M. (2018). On p-values and Bayes factors.
\emph{Annual Review of Statistics and Its Application}, \bold{5}, 393--419.

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
pCalibrate(p=c(0.05, 0.01, 0.001), type="exploratory")
zCalibrate(p=c(0.05, 0.01, 0.005), type="one.sided", 
           alternative="simple")
zCalibrate(p=c(0.05, 0.01, 0.005),  type="two.sided", 
           alternative="normal")
tCalibrate(p=c(0.05, 0.01, 0.005), n=c(10, 20, 50), type="two.sided", 
           alternative="normal")
FCalibrate(p=c(0.05, 0.01, 0.005), n=20, d=c(2, 5, 10), 
           alternative="chi.squared")
LRCalibrate(p=c(0.05, 0.01, 0.005), df=2, alternative="simple")
}
