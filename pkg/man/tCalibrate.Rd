\name{tCalibrate}
\alias{tCalibrate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Calibration of P-values from t-tests
%%  ~~function to do ... ~~
}
\description{
    Transforms one- and two-sided P-values from t-tests to sample-size adjusted lower bounds on the Bayes factor
    for the point null hypothesis against the alternative.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
tCalibrate(p, n, type="two.sided", alternative="normal", transform="id")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
    a vector of P-values
%%     ~~Describe \code{p} here~~
}
  \item{n}{
  a vector of sample sizes
%%     ~~Describe \code{n} here~~
}
  \item{type}{
   either \code{"one.sided"} or \code{"two.sided"}.
  Defaults to \code{"two.sided"}.
  Specifies if the t-test (and the corresponding P-value) is one-sided or two-sided.
}
  \item{alternative}{
  either \code{"simple"} or  \code{"normal"}.
  Defaults to \code{"normal"}.
  Specifies the alternative hypotheses for the mean to consider.
  See Details for more information.
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
  For one-sided t-tests (\code{type="one.sided"}), \code{alternative="simple"} considers all simple
 point alternatives for the mean and \code{alternative="normal"} all shifted and possibly mirrored (if \code{direction=less})
 half-normal distributions with mode at the null value.
 
  For two-sided t-tests (\code{type="two.sided"}), \code{alternative="simple"} considers all
  two-point distributions on the mean
  which are symmetric with respect to the null value and \code{alternative="normal"} 
  all local normal distributions centered around the null value.
  
  Note that for \code{alternative="normal"},
 the lower bound on the Bayes factor is the same for 
 one-sided and two-sided P-values.
 
    For \code{type="one.sided", alternative="simple"}, the calibration is available
    in closed form and is given in Held & Ott (2017), 
    equation (19).
    
    For \code{type="two.sided", alternative="simple"}, the calibration needs to be computed by
    numerical optimization and is given in Held & Ott (2017), 
    equations (17) and (18).
    
    For \code{alternative="normal"}, the corresponding lower bound on the Bayes factor has been derived in Johnson (2005) and the calibration is proposed in Held & Ott (2016), Section 3.
  
    If \code{alternative="normal"}, \code{tCalibrate()} is a special case of
    \code{FCalibrate()}, i.e.
     \code{tCalibrate(p, n, alternative="normal")}
     is equivalent to \cr
     \code{FCalibrate(p, n, d=1, alternative="normal")}.
    
    One can show that if one considers the class of all alternatives
which are symmetric with respect to the null value,
one obtains the same lower bound as for the subclass of all
symmetric two-point distributions.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
    A matrix containing the lower bounds on the Bayes factors 
    as entries,
    for all combinations of P-value and sample size.
    The values for the k-th sample size (k-th entry in the vector \code{n}) and 
    the specified P-values
    are given in the k-th row.
    %    (for the point null hypothesis against the composite alternative)
%    over all considered alternatives 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Held, L. and Ott, M. (2017). On p-values and Bayes factors.
 To appear in \emph{Annual Review of Statistics and Its Application, Volume 5}.

Held, L. and Ott, M. (2016). How the maximal evidence of P-values against point null hypotheses depends
on sample size. \emph{The American Statistician}, \bold{70}, 335--341.

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
    \code{\link{FCalibrate}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
tCalibrate(p=c(0.05, 0.01, 0.005), n=10, type="one.sided", alternative="simple")
# two-sided alternatives
tCalibrate(p=c(0.05, 0.01, 0.005), n=20, alternative="simple")
tCalibrate(p=c(0.05, 0.01, 0.005), n=c(10, 20, 50))

# for two-sided alternatives, plot the minimum Bayes factors 
# as a function of the P-value
par(mfrow=c(1,2), pty="s", las=1)

# plot for simple alternatives
p <- exp(seq(log(0.0001), log(0.3), by=0.01))
n <- c(10, 20, 50)
minBF <- tCalibrate(p, n, alternative = "simple")
# compare to the bound for large n
bound1 <- zCalibrate(p, alternative="simple")
matplot(p, t(minBF), ylim=c(0.0003, 1), type="l", ylab="Min. Bayes factor", log="xy", 
        xlab="two-sided t-test P-value", lty=1, lwd=2, axes=FALSE, 
        main="Simple alternative")
lines(p, bound1, col="gray", lty=2, lwd=2)
my.at1 <- c(0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3)
my.values1 <- as.character(c(format(c(0.0001,0.0003), nsmall=4, digits=4, scientific=FALSE), 
                    c(0.001, 0.003, 0.01, 0.03, 0.1, 0.3)))
axis(1, at=my.at1, labels=my.values1)
my.values2 <- c(3000, 1000, 300, 100, 30, 10, 3, 1)
my.at2 <- 1/my.values2
my.ylegend2 <- c(paste("1/", my.values2[-length(my.values2)], sep=""), "1")
axis(2, at=my.at2, labels=my.ylegend2)
box()
legend("bottomright", legend=rev(c("n=10", "n=20", "n=50", "n large")), 
        lty=rev(c(rep(1,3), 2)), lwd=2, col=rev(c(1:3, "gray")))
       
# plot for local normal alternatives
p <- exp(seq(log(0.0001), log(0.3), by=0.01))
n <- c(10, 20, 50)
minBF <- tCalibrate(p, n)
bound2 <- zCalibrate(p)
matplot(p, t(minBF), ylim=c(0.0003, 1), type="l", ylab="Min. Bayes factor", log="xy", 
        xlab="two-sided t-test P-value", lty=1, lwd=2, axes=FALSE, 
        main="Local normal alternatives")
lines(p, bound2, col="gray", lty=2, lwd=2)
axis(1, at=my.at1, labels=my.values1)
axis(2, at=my.at2, my.ylegend2)
box()            
legend("bottomright", legend=rev(c("n=10", "n=20", "n=50", "n large")), 
        lty=rev(c(rep(1,3), 2)), lwd=2, col=rev(c(1:3, "gray")))


}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
\keyword{models & regression}
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
