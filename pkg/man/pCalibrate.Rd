\name{pCalibrate}
\alias{pCalibrate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calibrations of two-sided p-values directly based on the p-value
%%  ~~function to do ... ~~
}
\description{
    Transforms a two-sided p-value to a minimum Bayes factor.
    That minimum Bayes factor is obtained by modelling the distributions
    of the p-value under the null and the alternative hypothesis,
    respectively.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pCalibrate(p,  type="exploratory", transform="id")
%pCalibrate(p,  alternative="noninformative", transform="id")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
  a vector of two-sided p-values
%%     ~~Describe \code{p} here~~
}
  \item{type}{
   either \code{"exploratory"} or \code{"confirmatory"}.
  Defaults to \code{"exploratory"}.
 % either \code{"noniformative"} or \code{"informative"}.
 % Defaults to \code{"noniformative"}.
  Corresponds to different distrbutional assumptions for the p-value p
  under the alternative.
}
    \item{transform}{
  either \code{"id"}, \code{"log"}, \code{"log2"}  or \code{"log10"}.
  Defaults to \code{"id"}.
  Specifies how to transform the minimum Bayes factor(s).
  \code{"id"} corresponds to no transformation.
  \code{"log"} refers to the natural logarithm, \code{"log2"} to the logarithm
  to the base 2 and \code{"log10"} to the logarithm
  to the base 10.
}
}
\details{
      If \code{type="exploratory"} is used,
      the so-called "- e p log(p)" calibration (Sellke et al., 2001)
      is applied.
      
      If \code{type="confirmatory"} is used,
      the  the so-called "- e q log(q)" calibration with q=1-p (Held & Ott, 2018)
      is applied.
  
      \code{type="exploratory"} gives a larger minimum Bayes factor than \code{type="confirmatory"}.
      
      Under the null hypothesis, the distribution of the p-value
      is assumed to be uniform.
      Under the alternative, the p-value is assumed to have a beta-distribution
      with monotonically decreasing density function under both alternatives,
      with different parameters however.
      If \code{type="exploratory"}, the prior sample size
      does not exceed 2, whereas for \code{type="confirmatory"},
      the prior sample size is at least 2.
      The latter calibration may be appropriate for small sample size,
      but for larger sample size it is too conservative (Held & Ott, 2018).
      
      Note that for the "- e p log(p)" calibration,
      alternative derivations which do not assume a beta-distribution under the
      alternative have also been given, see Sellke et al. (2001).
%%  ~~ If necessary, more details than the description above ~~
}
\note{
    The argument \code{type} replaces the argument \code{alternative}
    in version 0.1-1 of this package. \cr
    \code{type="exploratory"} corresponds to \code{alternative="noninformative"}
    and \cr
    \code{type="confirmatory"} corresponds to \code{alternative="informative"}.
}
\value{
A numeric vector of minimum Bayes factors for the specified p-values
}
\references{
Held, L. and Ott, M. (2018). On p-values and Bayes factors.
\emph{Annual Review of Statistics and Its Application}, \bold{5}, 393--419.
 
Sellke, T., Bayarri, M. J. and Berger, J. O. (2001). Calibration of p values for testing precise null hypotheses.
\emph{The American Statistician}, \bold{55}, 62--71.

Vovk, V. G. (1993). A logic of probability, with application to the foundations of statistics (with discussion and a reply by the author). \emph{Journal of the Royal Statistical Society, Series B}, \bold{55}, 317--351.
%% ~put references to the literature/web site here ~
}
%\author{
%Manuela Ott
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
pCalibrate(p=c(0.05, 0.01, 0.001))
pCalibrate(p=c(0.05, 0.01, 0.001), type="confirmatory")

# plot the 2 calibrations as a function of the p-value
par(las=1)
p <- exp(seq(log(0.0001), log(0.3), by=0.01))
minBF1 <- pCalibrate(p=p)
minBF2 <- pCalibrate(p=p, type="confirmatory")
plot(p, minBF1, type="l", log="xy",
      xlab="Two sided p-value", ylab="Minimum Bayes factor", 
     axes=FALSE, lwd=2, col=1)
lines(p, minBF2, col=2, lwd=2)
axis(1, at=c(0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3), 
     as.character(c(format(c(0.0001,0.0003), nsmall=4, digits=4, 
                  scientific=FALSE), 
     c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3))))
my.values <- c(3000, 1000, 300, 100, 30, 10, 3, 1)
my.at <- 1/my.values
my.ylegend <- c(paste("1/", my.values[-length(my.values)], sep=""), "1")
axis(2, at=my.at, my.ylegend)
box()
legend("bottomright", lty=1, lwd=2, 
        legend=c("- e p log(p)", "- e q log(q)")
                , col=c(1,2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
