\name{BF2pp}
\alias{BF2pp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Transform Bayes factors to posterior probabilities
%%  ~~function to do ... ~~
}
\description{
    Transforms the Bayes factor for a hypothesis H_1 against a hypothesis H_2
    to the posterior probability of H_1 given
    the prior probability of H_1 
    (assuming that the prior probabilities of H_1 and H_2 add up to 1).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BF2pp(BF, prior.prob=0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BF}{
   a vector of Bayes factors
%%     ~~Describe \code{bf} here~~
}
   \item{prior.prob}{
    a vector of prior probabilities for H_1.
    Defaults to a vector with entries 0.5.
%%     ~~Describe \code{prior.odds} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{A numeric vector of posterior probabilities of the hypothesis H_1
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
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
BF2pp(BF=c(2, 10, 50), prior.prob=c(0.2, 0.5, 0.9))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
