\name{FCalibrate}
\alias{FCalibrate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calibration of two-sided P-values from the F-test in the linear model
%%  ~~function to do ... ~~
}
\description{
      Transforms two-sided P-values from the F-test of overall significance in the linear model to sample-size adjusted lower bounds on the Bayes factor
    for the point null hypothesis against the alternative.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FCalibrate(p, n, d, alternative="chi.squared", intercept=TRUE, transform="id")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{a vector of two-sided P-values
%%     ~~Describe \code{p} here~~
}
  \item{n}{a scalar or a vector of positive integers. Specifies the sample size(s). May be a vector only if \code{d} is a scalar.
%%     ~~Describe \code{n} here~~
}
  \item{d}{a scalar or a vector of positive integers. Specifies the dimension(s) of the vector(s) of regression coefficients. May be a vector only if \code{n} is a scalar.
%%     ~~Describe \code{d} here~~
}
  \item{alternative}{either \code{"simple"} or \code{"chi.squared"}.
  Defaults to \code{"chi.squared"}.
  Specifies the alternative hypotheses
  on the non-centrality parameter of the F-distribution to consider.
  \code{"simple"} only considers simple point alternative hypotheses.
  \code{"chi.squared"} assumes a scaled chi-square distribution.
%%     ~~Describe \code{alternative} here~~
}
  \item{intercept}{logical. If TRUE, the linear model contains an unknown intercept term,
  otherwise the intercept is fixed. Defaults to TRUE.
}
    \item{transform}{
  either \code{"id"}, \code{"log"}, \code{"log2"}  or \code{"log10"}.
  Defaults to \code{"id"}.
  Specifies how to transform the lower bound on the Bayes factor.
  \code{"id"} corresponds to no transformation.
  \code{"log"} refers to the natural logarithm, \code{"log2"} to the logarithm
  to the base 2 and \code{"log10"} to the logarithm
  to the base 10.
}
}
\details{
  Note that under the point null hypothesis that all regression coefficients are
  equal to zero, the F-statistic \eqn{F}
  (which is the \eqn{(1-p)}-quantile of the F-distribution with \eqn{d} and \eqn{n-d-1} degrees of freedom)
  has a central F-distribution with \eqn{d} and \eqn{n-d-1} degrees of freedom
  if the linear model contains an unknown intercept term
  (otherwise \eqn{F}
  has a central F-distribution with \eqn{d} and \eqn{n-d} degrees of freedom).
  Under a simple point alternative, 
  F has a non-central F-distribution with \eqn{d} and \eqn{n-d-1} degrees of freedom. 
  
  To obtain the lower bound on the Bayes factor 
  for \code{alternative="simple"},
  the likelihood under the alternative is then maximized numerically with respect to
  the non-centrality parameter.
  That calibration is described in Held & Ott (2017), 
    Section 3.1 (in the last two paragraphs).

    
    For \code{alternative="chi-squared"}, the calibration is proposed in Held & Ott (2016), 
    Section 3
    and also described in Held & Ott (2017), Section 3.2.
    The corresponding minimum Bayes factor has already been derived in Johnson (2005).
    As described there, assigning a scaled chi-sqaured distribution to the
    non-centrality parameter of the F-distribution corresponds to
    assigning a (multivariate) normal prior distribution 
centered around the null value to the vector of regression coefficients.
}
\note{
    Computation may fail for \code{alternative="simple"} if the P-value \code{p} is extremely small
    and min\{n, n-d\} is also small. Warnings will be given in this case and
    the returned value is \code{minBF=NaN}.
}
\value{
     A matrix containing the lower bounds on the Bayes factors 
     as entries,
    for all combinations of P-value and sample size \code{n} or dimension \code{d}
    (whichever is multidimensional).
    The values for the k-th sample size or dimension (k-th entry in the vector \code{n} or \code{d}) and the different P-values
    are given in the k-th row.
}
\references{
Held , L. and Ott, M. (2016). How the maximal evidence of P -values against point null hypotheses depends on sample size. \emph{American Statistician}, \bold{70}, 335--341

Held, L. and Ott, M. (2017). On p-values and Bayes factors.
 To appear in \emph{Annual Review of Statistics and Its Application, Volume 5}.
 
Johnson, V. E. (2005). Bayes factors based on test statistics. 
 \emph{Journal of the Royal Statistical Society, Series B} \bold{67}:689--701.
%% ~put references to the literature/web site here ~
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{tCalibrate}} 
    % covers the special case d=1.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
FCalibrate(p=c(0.05, 0.01, 0.005), n=20, d=2, alternative="simple")
# chi-squared alternatives
FCalibrate(p=c(0.05, 0.01, 0.005), n=20, d=2, intercept=FALSE)
FCalibrate(p=c(0.05, 0.01, 0.005), n=20, d=c(2, 5, 10))
FCalibrate(p=c(0.05, 0.01, 0.005), n=c(10, 20, 50), d=2)

# plot for chi-squared alternatives: d=2 and different sample sizes n
# note that the minimum Bayes factor decreases with decreasing sample size
p <- exp(seq(log(0.0001), log(0.3), by=0.01))
n <- c(5, 10, 20)
minBF <- FCalibrate(p, n, d=2) 
# compare to the bound for large n
minTBF <- LRCalibrate(p, df=2)
par(las=1)
matplot(p, t(minBF), ylim=c(0.0003, 1), type="l", ylab="Minimum Bayes factor", log="xy", 
        xlab="two-sided F-test P-value", lty=1, lwd=2, axes=FALSE, 
        main="Local normal alternatives")
lines(p, minTBF, col="gray", lty=2, lwd=2)
axis(1, at=c(0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3), 
     as.character(c(format(c(0.0001,0.0003), nsmall=4, digits=4, scientific=FALSE), 
                    c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3))))
my.values <- c(3000, 1000, 300, 100, 30, 10, 3, 1)
my.at <- 1/my.values
my.ylegend <- c(paste("1/", my.values[-length(my.values)], sep=""), "1")
axis(2, at=my.at, my.ylegend)
box()
legend("bottomright", 
       legend=rev(c("n=5", "n=10", "n=20", "n large")), 
       col=rev(c(1:3, "gray")), lty=c(2, rep(1, times=3)), lwd=2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{models & regression}
% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
