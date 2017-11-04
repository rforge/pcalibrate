\name{LRCalibrate}
\alias{LRCalibrate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calibration of two-sided P-values obtained from the likelihood ratio test statistic
%%  ~~function to do ... ~~
}
\description{
    Transforms two-sided P-values from likelihood ratio (deviance) tests to lower bounds on the Bayes factor
    for the point null hypothesis against the alternative.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LRCalibrate(p, df, alternative="gamma", transform="id")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
   a vector of two-sided P-values
%%     ~~Describe \code{p} here~~
}
  \item{df}{
  a vector of degrees of freedom of the asymptotic chi-squared distribution(s)
  of likelihood ratio test statistic(s)
%%     ~~Describe \code{df} here~~
}
  \item{alternative}{
  either \code{"simple"} or \code{"gamma"}.
  Defaults to \code{"gamma"}.
  Specifies the alternative hypotheses on the non-centrality
  parameter of the non-central chi-squared distribution to consider.
  \code{"simple"} only considers simple point alternative hypotheses.
  \code{"gamma"} assumes a specific gamma distribution.
%%     ~~Describe \code{alternative} here~~
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
    Under the assumption that the parameter vector of interest
    (which has dimension \code{df}) is equal to the vector of zeros,
    the distribution of the deviance converges to a
    chi-squared distribution with \code{df} degrees of freedom.
    Under a simple point alternative for the parameter vector of interest
    and some regularity conditions,
     the distribution of the deviance converges to a
    non-central chi-squared distribution with \code{df} degrees of freedom.
    
    For \code{alternative = "simple"}, the lower bound on the Bayes factor
    is obtained by mazimizing the
    density function of the (asymptotic) non-central chi-squared distribution under the alternative
    with respect to the non-centrality parameter.
    That calibration is described in
    Held and Ott (2017), Section 4.2.1.
    
    The calibration for \code{alternative = "normal"} uses the test-based Bayes factors
    introduced in Johnson (2008).
    That approach is also outlined in Held and Ott (2017), Section 4.2.2.
    
    Using \code{alternative = "gamma"} yields a larger bound than
    \code{alternative = "simple"}.
    Typical applications of these calibrations include generalized linear models.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
    A matrix containing the lower bounds on the Bayes factors 
    as entries,
    for all combinations of P-value and degrees of freedom.
    The values for the k-th degrees of freedom (k-th entry in the vector \code{df}) and the different P-values
    are given in the k-th row.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Held, L. and Ott, M. (2017). On p-values and Bayes factors.
 To appear in \emph{Annual Review of Statistics and Its Application, Volume 5}.

Johnson, V. E. (2008). Properties of Bayes factors based on test statistics. \emph{Scandinavian
Journal of Statistics}, \bold{35}, 354--368.
%% ~put references to the literature/web site here ~
}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
LRCalibrate(p=c(0.05, 0.01, 0.005), df=2, alternative="simple")
# gamma alternatives
LRCalibrate(p=c(0.05, 0.01, 0.005), df=c(2, 5, 10))

# plot the minimum Bayes factor as a function of the P-value
# for different degrees of freedom df of the LR test statistic
par(mfrow=c(1,2), las=1)

p <- exp(seq(log(0.005), log(0.3), by=0.01))
df <- c(1, 5, 20)
par(las=1)

# for a simple alternative
minBF.sim <- LRCalibrate(p, df=df, alternative="simple")
matplot(p, t(minBF.sim), type="l", ylab="Minimum Bayes factor", log="xy", 
        xlab="Two-sided LR-test P-value", lty=1, lwd=2, axes=FALSE,
        main="Simple alternative")
axis(1, at=c(0.01, 0.03, 0.1, 0.3), c(0.01, 0.03, 0.1, 0.3))
my.values <- c(30, 20, 10, 5, 3, 1)
my.at <- 1/my.values
my.ylegend <- c(paste("1/", my.values[-length(my.values)], sep=""), "1")
axis(2, at=my.at, my.ylegend)
box()
legend("bottomright", legend=c("df=1", "df=5", "df=20"), 
        lty=1, lwd=2, , col=1:3)

# for gamma alternatives
minBF.loc <- LRCalibrate(p, df=df, alternative="gamma")
matplot(p, t(minBF.loc), type="l", ylab="  Minimum Bayes factor", log="xy", 
        xlab="Two-sided LR-test P-value", lty=1, lwd=2, axes=FALSE,
        main="Local alternatives")
axis(1, at=c(0.01, 0.03, 0.1, 0.3), c(0.01, 0.03, 0.1, 0.3))
axis(2, at=my.at, my.ylegend)
box()
legend("bottomright", legend=c("df=1", "df=5", "df=20"), 
        lty=1, lwd=2, col=1:3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{models & regression}
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
