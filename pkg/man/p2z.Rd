\name{p2z}
\alias{p2z}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Transform p-values to z-values
%%  ~~function to do ... ~~
}
\description{
    Transforms the p-value(s) to the corresponding one- or two-sided z-value(s).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
p2z(p, type="two.sided")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
   a vector of p-values
%%     ~~Describe \code{bf} here~~
}
   \item{type}{
    type of test. Either "two.sided" or "one.sided". 
    Defaults to "two.sided".
%%     ~~Describe \code{prior.odds} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{A numeric vector of z-values
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
p2z(p=c(0.005, 0.01, 0.05))
p2z(p=c(0.005, 0.01, 0.05), type="one.sided")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
