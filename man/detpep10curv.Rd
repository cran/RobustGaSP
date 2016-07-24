\name{detpep10curv}
\alias{detpep10curv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Dette & Pepelyshev (2010) Curved Function
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
A 3-dimensional test function.
%This function is used for the comparison of computer experiment designs. It is highly curved.
}
\usage{
detpep10curv(xx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xx}{ a 3-dimensional vector specifying the location where the function is to be evaluated.  

%%     ~~Describe \code{xx} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
See \url{http://www.sfu.ca/~ssurjano/detpep10curv.html} for more details.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A real value equal to the function value evaluated at xx.
}
\references{
%% ~put references to the literature/web site here ~
Dette, H., & Pepelyshev, A. (2010). Generalized Latin hypercube design for computer experiments. Technometrics, 52(4).

\url{http://www.sfu.ca/~ssurjano/detpep10curv.html}

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
%%##--  or do  help(data=index)  for the standard data sets.
%%
%%## The function is currently defined as
%%function (x) 
%%{
%%  }
%}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{internal}
