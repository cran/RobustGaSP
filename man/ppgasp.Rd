\name{ppgasp}
\alias{ppgasp}
%\alias{show.rgasp}
\alias{ppgasp-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Setting up the parallel partial GaSP model
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Setting up the parallel partial GaSP model for estimating the parameters (if the parameters are not given). 
%The range and noise-variance ratio parameters are given and/or have been estimated.
}
\usage{
   ppgasp(design,response,trend=matrix(1,dim(response)[1],1),zero.mean="No",nugget=0,
    nugget.est=F,range.par=NA,method='post_mode',prior_choice='ref_approx',a=0.2,
    b=1/(length(response))^{1/dim(as.matrix(design))[2]}*(a+dim(as.matrix(design))[2]),
    kernel_type='matern_5_2',isotropic=F,R0=NA,
    optimization='lbfgs',
    alpha=rep(1.9,dim(as.matrix(design))[2]),lower_bound=T,
    max_eval=max(30,20+5*dim(design)[2]),initial_values=NA,num_initial_values=2)
 % \S4method{show}{rgasp}(object)    
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{design}{ a matrix of inputs.
%%     ~~Describe \code{design} here~~
}
  \item{response}{ a matrix of outputs where each row is one runs of the computer model  output.
%%     ~~Describe \code{response} here~~
}
  \item{trend}{ the mean/trend matrix of inputs. The default value is a  vector of ones. 
%%     ~~Describe \code{trend} here~~
}
\item{zero.mean}{ it has zero mean or not. The default value is FALSE meaning the mean is not zero. TRUE means the mean is zero.
}
  \item{nugget}{ numerical value of the nugget variance ratio. If nugget is equal to 0, it means there is either no nugget or the nugget is estimated. If the nugget is not equal to 0, it means a fixed nugget. The default value is 0. 
%%     ~~Describe \code{nugget} here~~
}
  \item{nugget.est}{ boolean value. \code{T} means nugget should be estimated and \code{F} means nugget is fixed
  or not estimated. The default value is \code{F}.
%%     ~~Describe \code{nugget.est} here~~
}
  \item{range.par}{ either \code{NA} or a \code{vector}. If it is \code{NA}, it means range parameters are estimated; otherwise range parameters are given. The default value is \code{NA}.
%%     ~~Describe \code{range.par} here~~
}
\item{method}{  method of parameter estimation. \code{post_mode} means the marginal posterior mode is used for estimation. \code{mle} means the maximum likelihood estimation is used. \code{mmle} means the maximum marginal likelihood estimation is used. The \code{post_mode}  is the default method. 
}
  \item{prior_choice}{ the choice of prior for range parameters and noise-variance parameters. \code{ref_xi} and \code{ref_gamma} means the reference prior with reference prior with the log of inverse range parameterization \ifelse{html}{\out{&xi;}}{\eqn{\mathbf{\xi}}{xi}} or range parameterization \ifelse{html}{\out{&gamma;}}{\eqn{\mathbf{\gamma}}{gamma}}. \code{ref_approx} uses the jointly robust prior to approximate the reference prior. The default choice is \code{ref_approx}.
%%     ~~Describe \code{prior_choice} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
prior parameters in the jointly robust prior. The default value is 0.2. 
}
  \item{b}{
%%     ~~Describe \code{b} here~~
prior parameters in the jointly robust prior. The default value is \code{n^{-1/p}(a+p)} where n is the number of runs and p is the dimension of the input vector. 
}
  \item{kernel_type}{
%%     ~~Describe \code{kernel_type} here~~

A vector specifying the type of kernels of each coordinate of the input. \code{matern_3_2} and \code{matern_5_2} are \code{Matern correlation} with roughness parameter 3/2 and 5/2 respectively. \code{pow_exp} is power exponential correlation with roughness parameter alpha. If \code{pow_exp} is to be used, one needs to specify its roughness parameter alpha. The default choice is \code{matern_5_2}.
}

\item{isotropic}{
a boolean value. \code{T} means the isotropic kernel will be used and \code{F} means the separable kernel will be used. The default choice is  the separable kernel. 
}

\item{R0}{
the distance between inputs. If the value is \code{NA}, it will be calculated later. It can also be specified by the user. If specified by user, it is either a \code{matrix} or \code{list}. The default value is \code{NA}.
}
\item{optimization}{
the method for numerically optimization of the kernel parameters.  Currently three methods are implemented. \code{lbfgs} is the low-storage version of the Broyden-Fletcher-Goldfarb-Shanno method. \code{nelder-mead} is the  Nelder and Mead method.  \code{brent} is the Brent method for one-dimensional problems. 
 
}

  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
roughness parameters in the \code{pow_exp} correlation functions. The default choice is a vector with each entry being 1.9. 
}

 % \item{object}{An object of the corresponding signature.}  %%??
  \item{lower_bound}{
  boolean value.  \code{T} means the default lower bounds of the inverse range parameters are used to constrained the optimization and \code{F} means the optimization is unconstrained. The default value is \code{T} and we also suggest to use \code{F} in various scenarios. 
  }
    \item{max_eval}{
the maximum number of steps to estimate the range and nugget parameters.  
  }
  \item{initial_values}{
a matrix of initial values of the kernel parameters to be optimized numerically, where each row of the matrix contains a set of the log inverse range parameters and the log nugget parameter. 
  }
  \item{num_initial_values}{
the number of initial values of the kernel parameters in optimization.
  }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%expand here the details.
%}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\code{ppgasp} returns a S4 object of class \code{ppgasp} (see \code{ppgasp-class}).
}

\references{
%% ~put references to the literature/web site here ~
M. Gu. and J.O. Berger (2016). Parallel partial Gaussian process emulation for computer models with massive output. \emph{Annals of Applied Statistics}, 10(3), 1317-1347.

M. Gu, X. Wang and J.O. Berger (2018), \emph{Robust Gaussian stochastic process emulation}, \emph{Annals of Statistics}, 46(6A), 3038-3066.

M. Gu (2018), \emph{Jointly robust prior for Gaussian stochastic process in emulation, calibration and variable selection}, arXiv:1804.09329.

J. Nocedal (1980), Updating quasi-Newton matrices with limited storage, \emph{Math. Comput.}, 35, 773-782.


D. C. Liu and J. Nocedal (1989), On the limited memory BFGS method for large scale optimization, \emph{Math. Programming}, 45, p. 503-528.

Brent, R. (1973), Algorithms for Minimization without Derivatives. Englewood Cliffs N.J.: Prentice-Hall.


}

\author{
\packageAuthor{RobustGaSP}

Maintainer: \packageMaintainer{RobustGaSP}
}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{

  library(RobustGaSP)
  
  ###parallel partial Gaussian stochastic process (PP GaSP) model 
  ##for the humanity model
  data(humanity_model)
  ##120 runs. The input has 13 variables and output is 5 dimensional.
  ##PP GaSP Emulator
  m.ppgasp=ppgasp(design=humanity.X,response=humanity.Y,nugget.est= TRUE)
  show(m.ppgasp)

  ##make predictions
  m_pred=predict(m.ppgasp,humanity.Xt)
  sqrt(mean((m_pred$mean-humanity.Yt)^2))
  mean(m_pred$upper95>humanity.Yt & humanity.Yt>m_pred$lower95)
  mean(m_pred$upper95-m_pred$lower95)
  sqrt( mean( (mean(humanity.Y)-humanity.Yt)^2 ))

  ##with a linear trend on the selected input performs better
  \dontrun{
    ###PP GaSP Emulation with a linear trend for the humanity model
    data(humanity_model)
    ##pp gasp with trend
    n<-dim(humanity.Y)[1]
    n_testing=dim(humanity.Yt)[1]
    H=cbind(matrix(1,n,1),humanity.X$foodC)
    H_testing=cbind(matrix(1,n_testing,1),humanity.Xt$foodC)
    m.ppgasp_trend=ppgasp(design=humanity.X,response=humanity.Y,trend=H, 
    nugget.est= TRUE)
    
    show(m.ppgasp_trend)
    
    ##make predictions
    m_pred_trend=predict(m.ppgasp_trend,humanity.Xt,testing_trend=H_testing)
    sqrt(mean((m_pred_trend$mean-humanity.Yt)^2))
    mean(m_pred_trend$upper95>humanity.Yt & humanity.Yt>m_pred_trend$lower95)
    mean(m_pred_trend$upper95-m_pred_trend$lower95)
  }

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
