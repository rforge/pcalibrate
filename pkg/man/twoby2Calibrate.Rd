\name{twoby2Calibrate}
\alias{twoby2Calibrate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Minimum Bayes factors and P-values from Fisher's exact test for 2x2 contingency tables
%%  ~~function to do ... ~~
}
\description{
    Computes a sample-size adjusted minimum Bayes factor
    for the given 2x2 contingency table.
    Also returns P-values from Fisher's exact test
    (different versions in the two-sided case) and less conservative
    alternatives such as a mid P-value (see Details for more information).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
twoby2Calibrate(x, type="two.sided", alternative="normal", direction=NULL, 
                transform.bf="id")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a 2x2 contingency table in matrix form
%%     ~~Describe \code{tab} here~~
}
   \item{type}{
   either \code{"one.sided"} or \code{"two.sided"}.
  Defaults to \code{"two.sided"}.
  Specifies if Fisher's exact test (and the corresponding P-value) is one-sided or two-sided.
}
  \item{alternative}{
   either \code{"simple"} or \code{"normal"}.
  Defaults to \code{"normal"}.
  Specifies the alternative hypotheses for the (log) odds ratio to consider
  for two-sided tests.
  Is ignored if \code{type="one.sided"}
  (in this case only simple alternatives are available).
%%     ~~Describe \code{alternative} here~~
}
\item{direction}{
   either \code{"greater"}, \code{"less"} or \code{NULL}.
   Defaults to \code{NULL}.
  Specifies the direction of the alternative for one-sided tests: \code{"greater"} corresponds
  to an odds ratio > 1 and \code{"less"} to an odds ratio < 1.
  Is ignored if \code{type="two.sided"}.
}
    \item{transform.bf}{
  either \code{"id"}, \code{"log"}, \code{"log2"}  or \code{"log10"}.
  Defaults to \code{"id"}.
  Specifies how to transform the minimum Bayes factor.
  \code{"id"} corresponds to no transformation.
  \code{"log"} refers to the natural logarithm, \code{"log2"} to the logarithm
  to the base 2 and \code{"log10"} to the logarithm
  to the base 10.
}
}
\details{
  If \code{type="two.sided"},
  the point null hypothesis that the odds ratio is 1 is tested
  against
  specific two-sided alternatives:
  \code{alternative="simple"} considers all two-point distributions symmetric around 0 for
  the log odds ratio.
  \code{alternative="normal"} assumes a local normal prior distribution (a so-called g-prior) centered around 0
  for the log odds ratio.  
  
  In the one-sided case (\code{type="one.sided"}), \code{direction="less"} tests the alternative that the odds ratio is
  less than 1 and considers simple point alternatives in that direction
  to compute the minimum Bayes factor.
  \code{"greater"} does the same for the alternative that the odds ratio is
  larger than 1.
  
  The calibration obtained with \code{type="two.sided", alternative="normal"} is based on
  the methodology proposed in Li & Clyde (2016) 
  and yields an (approximate) minimum Bayes factor
  in closed form.
  All the other minimum Bayes factors
  are computed by numerical optimization.
  For \code{type="two.sided"},
  the two calibrations are described in Ott & Held (2017).
  
  For one-sided alternatives, the following P-value and 2 related quantities are computed:
  \itemize{
   \item{ \code{p.fi} }{  is the one-sided P-value from Fisher's exact test.}
   \item{ \code{p.mid} }{is a "mid" P-value. It is obtained by subtracting half of the
   probability mass of the observed table from \code{p.fi}.}
   \item{ \code{p.lie} }{is a Bayesian posterior probability. 
   If \code{direction="greater"}, it is the posterior probability that the odds ratio
   exceeds 1 given the observed table under the assumption of uniform priors on the
   success probabilities for the two groups.
   If \code{direction="less"}, it is the posterior probability that the odds ratio
   does not exceed 1 given the observed table under the same priors.}
   }
   
  For two-sided alternatives, the following 3 P-values and 2 related quantities are computed:
  \itemize{
   \item{ \code{p.pb} }{  is the "probability-based" P-value (the classical choice), 
   defined as the sum
   of the probabilities of all tables which are at most a likely as the observed table and
   have the same marginals.}
   \item{ \code{p.ce} }{ is the "central" P-value, which is twice the minimum one-sided P-value
    (from Fisher's exact test),
   bounded by 1.}
   \item{ \code{p.bl} }{is "Blaker's" P-value, which is the minimum one-sided P-value
   (from Fisher's exact test)
   plus the largest tail probability from the other tail of the distribution that
   does not exceed that minimum.}
   \item{ \code{p.mid} }{is a "mid" P-value. It is the mid-p modification of the central P-value, 
   i.e. it equals twice the minimum
   one-sided mid P-value.} 
   \item{ \code{p.lie} }{is a two-sided version of the posterior probability
    for the one-sided test. 
   Let p.lie.os be the one-sided posterior probability that the odds ratio does not exceed 1
   given the observed table,
   as returned by the one-sided test with \code{direction="less"}.
   Then \code{p.lie} = 2 min\{p.lie.os, 1-p.lie.os\}.}
   }
   
   For one-sided alternatives,
   the posterior probability \code{p.lie} was already studied in Liebermeister (1877)
 and its frequentist properties are investigated in Seneta & Phipps (2001).
   
  For two-sided alternatives, the 3 P-values from Fisher's exact test are defined in equations (2.24)-(2.26)
in Kateri (2014) and computed using the \code{exact2x2()}-function in the 
\code{exact2x2}-package.
  The "mid" P-value is described in Rothman & Greenland (1998, pp. 222-223)
and computed using the \code{tab2by2.test()}-function
in the package \code{epitools}.
 The Bayesian significance measure \code{p.lie} is proposed in Ott & Held (2017)
 as a modification of the corresponding one-sided significance measure.
}
\value{A list of the following two elements:
\item{minBF }{the minimum Bayes factor}
\item{p.value }{A vector of 3 one-sided P-values/significance measures
for one-sided tests,
namely first the P-value \code{p.fi} from Fisher's exact test,
second the corresponding mid P-value \code{p.mid} and 
 third the Bayesian posterior probability \code{p.lie}
 (see Details for more information).
 
A vector of 5 two-sided P-values/significance measures
for two-sided tests: 
The first three P-values
\code{p.pb}, \code{p.ce} and \code{p.bl} (see Details for the definitions) correspond to
two-sided P-values from Fisher's exact test.
The 4th quantity \code{p.mid} is a mid P-value, namely the mid-p modification of the second P-value \code{p.ce}.
The last element \code{p.lie} is a Bayesian significance measure
(see Details for additional information).}
}
\section{Warning}{
  For 2x2 tables with entries equal to 0,
  the minimum Bayes factor is either not defined 
  (for \code{alternative="normal"}) or 
  the underlying numerical optimization is unstable 
  (for \code{alternative="simple"}).
  A warning is be displayed in such cases and \code{minBF=NA}
  is returned,
  but the different P-values/significance measures 
  are still available.
}
\references{
Li, Y. and Clyde, M. A. (2016). Mixtures of g-priors in generalized linear models. Technical report, Clemson/Duke
University. \url{https://arxiv.org/abs/1503.06913v2}

Liebermeister, C. (1877). Ueber Wahrscheinlichkeitsrechnung in Anwendung auf therapeutische Statistik. Sammlung klinischer Vortraege. \emph{Innere Medizin}, \bold{110}(31--64), 935--
962.

Kateri, M. (2014). Contingency Table Analysis - Methods and Implementation using R.
Statistics for Industry and Technology. Birkhaeuser.

Ott, M. and Held, L. (2017). Bayesian calibration of P-values from Fisher's
exact test. Technical report, University of Zurich.

Rothman, K. J. and Greenland, S. (1998). Modern Epidemiology. 2nd ed. Lippincott-Raven.

Seneta, E. and Phipps, M. C. (2001). On the comparison of two observed frequencies.
\emph{Biometrical Journal}, \bold{43}(1), 23--43.
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
  For computation of P-values:
  \code{\link[exact2x2]{exact2x2}} in package exact2x2 and
  \code{\link[epitools]{tab2by2.test}} in package epitools.
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
tab <- matrix(c(1,15,5,10), nrow=2, byrow=TRUE)
minBF.plus <- twoby2Calibrate(x=tab, type="one.sided", direction="greater")$minBF
minBF.minus <- twoby2Calibrate(x=tab, type="one.sided", direction="less")$minBF
minBF.sim <- twoby2Calibrate(x=tab, type="two.sided", alternative="simple")$minBF
minBF.nor <- twoby2Calibrate(x=tab)$minBF
p.plus <- twoby2Calibrate(x=tab, type="one.sided", direction="greater")$p.value
p.minus <- twoby2Calibrate(x=tab, type="one.sided", direction="less")$p.value
pvals.twosid <- twoby2Calibrate(x=tab)$p.value
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\keyword{ htest }
\keyword{models & regression}
