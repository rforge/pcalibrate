\name{formatBF}
\alias{formatBF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Formatting of Bayes factors
%%  ~~function to do ... ~~
}
\description{
      Bayes factors <1 are represented as ratios
     1/x, where x is rounded to the specified number of digits.
     Bayes factors >=1 are rounded to the specified number of digits.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
formatBF(BF, digits="default")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{BF}{ a numeric vector of Bayes factors
%%     ~~Describe \code{BF} here~~
}
  \item{digits}{ either "default" (see Details) or a positive integer specifiying the number of decimal places to round the Bayes factor (for Bayes factors >= 1) or
     its inverse (for Bayes factors < 1) to
%%     ~~Describe \code{BF} here~~
}
}
\details{
      The default formatting, which  is recommended in Held and Ott (2017), is as follows :
      For very small Bayes factors BF < 1/1000, "< 1/1000" is returned.
    Bayes factors BF with 1/1000 <= BF <=1/10 are formatted as 1/x, where x is an integer
    and Bayes factors BF with 1/10 < BF < 1 as
    1/x, where x is rounded to one decimal place.
    Accordingly, Bayes factors 1 <= BF < 10 are rounded to one decimal place,
     Bayes factors 10 <= BF <= 1000 are rounded to the next integer
     and for larger Bayes factors, "> 1000" is returned.
     
     If \code{digits} is specified,
     the Bayes factor (if it is >= 1) or its inverse (if the Bayes factor is <1)
     is rounded to the number of decimal places specified
     and returned as a ratio if the Bayes factor is <1.
%%  ~~ If necessary, more details than the description above ~~
}
\value{  A character vector of ratios (for inputs <1) or rounded numeric values (for inputs >=1).
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
 Held, L. and Ott, M. (2017). On p-values and Bayes factors.
 To appear in \emph{Annual Review of Statistics and Its Application, Volume 5}.
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
 BF <- c(0.0001, 0.0049, 0.0258, 0.24, 2.798)
 # use the default formatting
 formatBF(BF)
 # specify the number of digits
 formatBF(c(0.087, 4.65), digits=1)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
