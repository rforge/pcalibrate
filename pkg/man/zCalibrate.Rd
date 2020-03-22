\name{zCalibrate}
\alias{zCalibrate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Calibration of p-values from z-tests
%%  ~~function to do ... ~~
}
\description{
    Transforms one- and two-sided p-values from z-tests to minimum Bayes factors.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
zCalibrate(p=NULL, z=p2z(p), type="two.sided",  
           alternative="normal", transform="id")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{a vector of p-values, default is \code{NULL}
%%     ~~Describe \code{p} here~~
}
  \item{z}{a vector of z-values, is calculated from \code{p}, if necessary.
%%     ~~Describe \code{p} here~~
}
   \item{type}{
   either \code{"one.sided"} or \code{"two.sided"}.
  Defaults to \code{"two.sided"}.
  Specifies if the z-test (and the corresponding p-value) is one-sided or two-sided.
}
  \item{alternative}{either \code{"simple"}, \code{"normal"} or \code{"local"}.
  Defaults to \code{"normal"}.
  Specifies the alternative hypotheses for the mean to consider.
  See Details for more information.
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
 For one-sided z-tests (\code{type="one.sided"}), 
 \code{alternative="simple"} considers all simple
 point alternatives for the mean in the specified direction, \code{alternative="normal"} all shifted and possibly mirrored (if \code{direction=less})
 half-normal distributions with mode at the null value and 
 \code{alternative="local"} all monotonic density functions with mode at the null value,
 which are concentrated on the parameter space of the alternative. 
 
  For two-sided z-tests (\code{type="two.sided"}), \code{alternative="simple"} considers all
  two-point distributions on the mean
  which are symmetric with respect to the null value, \cr \code{alternative="normal"} 
  all normal distributions centered around the null value and 
 \cr \code{alternative="local"} all unimodal symmetric distributions with respect the null value
 (i.e. these distributions are non-increasing as a function of the distance to the null value).
 

 Note that for \code{alternative="normal"} and \code{alternative="local"},
 the minimum Bayes factor is the same for \code{type="one.sided"} and 
 \code{type="two.sided"}.
 
For \code{type="one.sided", alternative="simple"} and for \code{alternative="normal"},
the calibrations have closed-form expressions and were proposed in Edwards et al. (1963).\cr
  For \code{type="two.sided", alternative="simple"}  and for \code{alternative="local"}, the calibrations need to be computed by numerical optimization and were derived in 
Berger & Sellke (1987). 
Most of these calibrations are also described in Held & Ott (2018).


Note that if one considers the class of all alternatives,
the resulting minimum Bayes factor is the same as for 
\code{type="one.sided", alternative="simple"}.

One can show that if one considers the class of all alternatives
which are symmetric around the null value,
one obtains the same minimum Bayes factor as for the subclass of all
symmetric two-point distributions (Berger & Sellke, 1987).
}
\value{
    A numeric vector of minimum Bayes factors for the specified p-values 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Berger, J. O. and Sellke, T. (1987). Testing a point null hypothesis: The irreconcilability of P
values and evidence (with discussion). \emph{Journal of the American Statistical Association}, \bold{82}, 112--139.

    Edwards, W., Lindman, H. and Savage, L. J. (1963). Bayesian statistical inference for psychological research. 
    \emph{Psychological Review}, \bold{70}, 193--242.
    
Held, L. and Ott, M. (2018). On p-values and Bayes factors.
\emph{Annual Review of Statistics and Its Application}, \bold{5}, 393--419.
%% ~put references to the literature/web site here ~
}
%\author{
%Manuela Ott
%}
\note{
For small p-values, using the corresponding z-value
as an argument instead of the p-value may lead to more accurate computations.
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{tCalibrate}} 
 for sample-size adjusted calibrations of p-values from t-tests
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{  
# two-sided alternatives
zCalibrate(p=c(0.05, 0.01, 0.005), alternative="simple")
zCalibrate(p=c(0.05, 0.01, 0.005))
zCalibrate(p=c(0.05, 0.01, 0.005), alternative="local")
zCalibrate(z=c(2, 3, 4), alternative="local")
# one-sided alternatives
zCalibrate(p=c(0.05, 0.01, 0.005), type="one.sided", 
           alternative="simple")
zCalibrate(p=c(0.05, 0.01, 0.005), type="one.sided")
zCalibrate(p=c(0.05, 0.01, 0.005), type="one.sided", 
           alternative="local")

# plot the different calibrations as a function of the p-value
par(las=1)
p <- exp(seq(log(0.0001), log(0.3), by=0.01))
minBF <- matrix(NA, ncol=4, nrow=length(p))
minBF[ ,1] <- zCalibrate(p)
minBF[ ,2] <- zCalibrate(p, type="two.sided", alternative="local")
minBF[ ,3] <- zCalibrate(p, type="two.sided", alternative="simple")
minBF[ ,4] <- zCalibrate(p, type="one.sided", alternative="simple")
matplot(p, minBF, type="l", ylab="Minimum Bayes factor", log="xy", 
        xlab="z-test p-value", lty=1, lwd=2, col=c(1,3,2,4), axes=FALSE)
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
        legend=c("normal", "unimodal symmetric", 
                "two-sided simple", "one-sided simple"), col=c(1,3,2,4))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
