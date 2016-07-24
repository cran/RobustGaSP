\name{Borehole}
\alias{borehole}
\alias{Borehole}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Borehole function
%%  ~~function to do ... ~~
}
\description{

It is an 8-dimensional test function that models water flow through a borehole. Its simplicity and quick evaluation makes it a commonly used function for testing a wide variety of methods in computer experiments. 

}
\usage{
borehole(xx)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xx}{ an 8-dimensional vector specifying the location where the function is to be evaluated. }
%%     ~~Describe \code{x} here~~
}

\details{
%%  ~~ If necessary, more details than the description above ~~
See \url{http://www.sfu.ca/~ssurjano/borehole.html} for more details.
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The response is water flow rate, in m^3/yr, evaluated at \code{xx}.


}
\references{
%% ~put references to the literature/web site here ~

Harper, W. V., & Gupta, S. K. (1983). Sensitivity/uncertainty analysis of a borehole scenario comparing Latin Hypercube Sampling and deterministic sensitivity approaches (No. BMI/ONWI-516). Battelle Memorial Inst., Columbus, OH (USA). Office of Nuclear Waste Isolation.

Worley, B. A. (1987). Deterministic uncertainty analysis (No. CONF-871101-30). Oak Ridge National Lab., TN (USA).

\url{http://www.sfu.ca/~ssurjano/borehole.html}


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
