\name{rgasp}
\alias{rgasp}
%\alias{show.rgasp}
\alias{rgasp-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Setting up the robust GaSP model
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Setting up the robust GaSP model for estimating the parameters (if the parameters are not given). 
%The range and noise-variance ratio parameters are given and/or have been estimated.
}
\usage{
  rgasp(design, response, trend = matrix(1, length(response), 1), 
      nugget = 0, nugget.est = F, range.par = NA, prior_choice, 
      a = 0.2, b = 1/(length(response))^{1/dim(as.matrix(design))[2]}*(a+
      dim(as.matrix(design))[2]), 
      kernel_type = "matern_5_2", 
      alpha = rep(1.9, dim(as.matrix(design))[2]),
      multiple_starts=T,lower_bound=T,max_eval=30,xtol_rel=1e-5)
 % \S4method{show}{rgasp}(object)    
      
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{design}{ a matrix of inputs.
%%     ~~Describe \code{design} here~~
}
  \item{response}{ a matrix of outputs.
%%     ~~Describe \code{response} here~~
}
  \item{trend}{ the mean/trend matrix of inputs. The default value is a  vector of ones. 
%%     ~~Describe \code{trend} here~~
}
  \item{nugget}{ numerical value of the nugget variance ratio. If nugget is equal to 0, it means there is either no nugget or the nugget is estimated. If the nugget is not equal to 0, it means a fixed nugget. The default value is 0. 
%%     ~~Describe \code{nugget} here~~
}
  \item{nugget.est}{ boolean value. \code{T} means nugget should be estimated and \code{F} means nugget is fixed
  or not estimated. The default value is F \code{F}.
%%     ~~Describe \code{nugget.est} here~~
}
  \item{range.par}{ either \code{NA} or a \code{vector}. If it is \code{NA}, it means range parameters are estimated; otherwise range parameters are given. The default value is \code{NA}.
%%     ~~Describe \code{range.par} here~~
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
type of kernel. \code{matern_3_2} and \code{matern_5_2} are \code{Matern correlation} with roughness parameter 3/2 and 5/2 respectively. \code{pow_exp} is power exponential correlation with roughness parameter alpha. If \code{pow_exp} is to be used, one needs to specify its roughness parameter alpha. The default choice is \code{matern_5_2}.
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
roughness parameters in the \code{pow_exp} correlation functions. The default choice is a vector with each entry being 1.9. 
}
 % \item{object}{An object of the corresponding signature.}  %%??
  \item{multiple_starts}{
  boolean value. \code{T} means it uses multiple default initial points to search the optimal range and nugget parameters  and \code{F} means it use one default points to search the optimal parameters. The default value is \code{T}. 
  }
  \item{lower_bound}{
  boolean value.  \code{T} means the default lower bounds of the inverse range parameters are used to constrained the optimization and \code{F} means the optimization is unconstrained. The default value is \code{T} and we also suggest to use \code{F} in various scenarios. 
  }
  \item{max_eval}{
  the maximum number of steps to estimate the range and nugget parameters. The default choice is 30. 
  }
  \item{xtol_rel}{
  the torerance bound to stop the searching algorithm in estimation of the range and nugget parameters. The default choice is \code{10^{-5}}.
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
\code{rgasp} returns a S4 object of class \code{rgasp} (see \code{rgasp-class}).
}

\references{
%% ~put references to the literature/web site here ~
Mengyang Gu. (2016). Robust Uncertainty Quantification and Scalable Computation for Computer Models with Massive Output. Ph.D. thesis. Duke University.

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
%\dontrun{
library(RobustGaSP)
 #------------------------
  # a 3 dimensional example
  #------------------------
  # dimensional of the inputs
  dim_inputs <- 3    
  # number of the inputs
  num_obs <- 30       
  # uniform samples of design
  input <- matrix(runif(num_obs*dim_inputs), num_obs,dim_inputs) 
  
  # Following codes use maximin Latin Hypercube Design, which is typically better than uniform
  # library(lhs)
  # input <- maximinLHS(n=num_obs, k=dim_inputs)  ##maximin lhd sample
  
  ####
  # outputs from the 3 dim detpep10curv function
  
  output = matrix(0,num_obs,1)
  for(i in 1:num_obs){
    output[i]<-detpep10curv (input[i,])
  }
  
  # use constant mean basis, with no constraint on optimization
  m1<- rgasp(design = input, response = output, lower_bound=FALSE)
  
  # the following use constraints on optimization
  # m1<- rgasp(design = input, response = output, lower_bound=TRUE)
  
  # the following use a single start on optimization
  # m1<- rgasp(design = input, response = output, lower_bound=FALSE, multiple_starts=FALSE)

show(m1)      ####show this rgasp object 

m1@beta       ####estimated inverse range parameters

  #-----------------------------------
  # a 2 dimensional example with trend
  #-----------------------------------
  # dimensional of the inputs
  dim_inputs <- 2    
  # number of the inputs
  num_obs <- 20       
  
  # uniform samples of design
  input <-matrix(runif(num_obs*dim_inputs), num_obs,dim_inputs) 
  # Following codes use maximin Latin Hypercube Design, which is typically better than uniform
  # library(lhs)
  # input <- maximinLHS(n=num_obs, k=dim_inputs)  # maximin lhd sample
  
  # outputs from the 2 dim Brainin function
  
  output <- matrix(0,num_obs,1)
  for(i in 1:num_obs){
    output[i]<-limetal02non (input[i,])
  }
  
  ####trend or mean basis
  X<-cbind(rep(1,num_obs), input )
  
  
  # use constant mean basis with trend, with no constraint on optimization
  m2<- rgasp(design = input, response = output,trend =X,  lower_bound=FALSE,multiple_start=TRUE)

show(m1)      # show this rgasp object 

m1@beta       # estimated inverse range parameters
m1@theta_hat  

  #--------------------------------------------------------------------------------------
  # an 8 dimensional example using only a subset inputs and a noise with unknown variance
  #--------------------------------------------------------------------------------------
  # dimensional of the inputs
  dim_inputs <- 8    
  # number of the inputs
  num_obs <- 30       
  
  # uniform samples of design
  input <-matrix(runif(num_obs*dim_inputs), num_obs,dim_inputs) 
  # Following codes use maximin Latin Hypercube Design, which is typically better than uniform
  # library(lhs)
  # input <- maximinLHS(n=num_obs, k=dim_inputs)  # maximin lhd sample
  
  # rescale the design to the domain
  input[,1]<-0.05+(0.15-0.05)*input[,1];
  input[,2]<-100+(50000-100)*input[,2];
  input[,3]<-63070+(115600-63070)*input[,3];
  input[,4]<-990+(1110-990)*input[,4];
  input[,5]<-63.1+(116-63.1)*input[,5];
  input[,6]<-700+(820-700)*input[,6];
  input[,7]<-1120+(1680-1120)*input[,7];
  input[,8]<-9855+(12045-9855)*input[,8];
  
  # outputs from the 8 dim Borehole function
  
  output=matrix(0,num_obs,1)
  for(i in 1:num_obs){
    output[i]=borehole(input[i,])
  }
  
  
    
    
  
  # use constant mean basis with trend, with no constraint on optimization
  m3<- rgasp(design = input[,c(1,4,6,7,8)], response = output, 
            nugget.est=TRUE, lower_bound=FALSE,multiple_start=TRUE)

m3@beta       # estimated inverse range parameters
m3@nugget     


%}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
