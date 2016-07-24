\name{limetal02non}
\alias{limetal02non}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Lim et at. (2002) nonpolynomial function
%%  ~~function to do ... ~~
}
\description{

This function is 2-dimensional test function. 
%It is an example of a nonpolynomial model which exhibits a shape similar to that of a multivariate polynomial. 

}
\usage{
limetal02non(xx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xx}{ a 2-dimensional vector specifying the location where the function is to be evaluated.}
%%     ~~Describe \code{x} here~~
}

\details{
See \url{http://www.sfu.ca/~ssurjano/limetal02non.html} for more details.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A real value of this function evaluated at \code{xx}.

}
\references{
%% ~put references to the literature/web site here ~
\url{http://www.sfu.ca/~ssurjano/limetal02non.html}

Lim, Y. B., Sacks, J., Studden, W. J., & Welch, W. J. (2002). Design and analysis of computer experiments when the output is highly correlated over the input space. , \emph{Canadian Journal of Statistics}, \bold{30}(1), 109-126.

}
\author{
%%  ~~who you are~~
Copyright 2013. Derek Bingham, Simon Fraser University.
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

%\examples{
%%##---- Should be DIRECTLY executable !! ----
%%##-- ==>  Define data, use random,
%%##--	or do  help(data=index)  for the standard data sets.
%%
%%## The function is currently defined as
%%function (x) 
%%{
%%  }
%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{optimize}
\keyword{internal}% __ONLY ONE__ keyword per line
