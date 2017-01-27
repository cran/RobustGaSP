##########################################################################
## Functions needed for the Robust GaSP
## 
## Robust GaSP Package
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 3, April 2013.
##
## Copyright (C) 2015-present Mengyang Gu, James O. Berger, Jesus Palomo 
##							  
##    
##########################################################################

####################
findInertInputs<-function(object,threshold=0.1){
  P_hat=object@p*object@beta_hat*object@CL/sum(object@beta_hat*object@CL)
  index_inert=which(P_hat<threshold)
  
  cat('The estimated normalized inverse range parameters are :', P_hat,'\n')
  if(length(which(P_hat<0.1))>0){
    cat('The inputs ', index_inert, 'are suspected to be inert inputs','\n')
  }else{
    cat('no input is suspected to be an inert input', '\n')
  }
  P_hat
}

####################
neg_log_marginal_post_approx_ref <- function(param,nugget, nugget.est,R0,X,output,CL,a,b,kernel_type,alpha) {
   #####this has mean X, we should also include the case where X is not zero
   #####
   lml=log_marginal_lik(param,nugget,nugget.est,R0,X,output,kernel_type,alpha);
   lp=log_approx_ref_prior(param,nugget,nugget.est,CL,a,b);
     
  -(lml+lp)

}

####################
neg_log_marginal_post_ref<- function(param,nugget, nugget.est,R0,X,output,prior_choice,kernel_type,alpha) {
  
  lmp=log_ref_marginal_post(param,nugget,nugget.est,R0,X,output,kernel_type,alpha);
  
  
  if(prior_choice=='ref_xi'){
     if(nugget.est==T){###note that in this case nu is also needed to be transformed have tail properties
      #-sum(param)-lmp  ###this will let nu be non zero so we need to be careful
       -sum(param[1:(length(param))])-lmp ###this might be fine when intercept is in the trend
     }else{ ###no nugget
      -sum(param)-lmp
     }
  }else if(prior_choice=='ref_gamma'){
    if(nugget.est==T){###note that in this case nu is also needed to be transformed have tail properties
      -2*sum(param[1:(length(param)-1)])-lmp 
    }else{ ###no nugget
      -2*sum(param)-lmp
    }
  }
}



########################
neg_log_marginal_post_approx_ref_deriv<- function(param,nugget,nugget.est,R0,X,output,CL,a,b,kernel_type,alpha) {
  lml_dev=log_marginal_lik_deriv(param,nugget,nugget.est,R0,X,output,kernel_type,alpha)
  lp_dev=log_approx_ref_prior_deriv(param,nugget,nugget.est,CL,a,b)
  
  -(lml_dev+lp_dev)*exp(param)
    
}
#############this is a function to search the lower bounds for range parameters beta
#####need R0 in the function
search_LB_prob<-function(param, R0,COND_NUM_UB,p,kernel_type,alpha,nugget){
  num_obs=dim(R0[[1]])[1]
  propose_prob=exp(param)/(exp(param)+1)
  LB=NULL
  for( i_LB in 1:p){
    LB=c(LB, log(-log(propose_prob)/(max(R0[[i_LB]]))))    ###LB is log beta
  }
  
  R=separable_kernel(R0,exp(LB),kernel_type,alpha)  
  
  R=as.matrix(R)
  R=R+nugget*diag(num_obs)
  (kappa(R)-COND_NUM_UB)^2
}
###############a 2 dim hig02 function
hig02 <- function(s)
{
  ##########################################################################
  #
  # HIGDON (2002) FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  #########################################################################
  
  term1 <- sin(2*pi*s/10)
  term2 <- 0.2 * sin(2*pi*s/2.5)
  
  y <- term1 + term2
  return(y)
}


#############a two dimension lim function
limetal02non <- function(xx)
{
  ##########################################################################
  #
  # LIM ET AL. (2002) FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUT:
  #
  # xx = c(x1, x2)
  #
  ##########################################################################
  
  x1 <- xx[1]
  x2 <- xx[2]
  
  fact1 <- 30 + 5*x1*sin(5*x1)
  fact2 <- 4 + exp(-5*x2)
  
  y <- (fact1*fact2 - 100) / 6
  return(y)
}

#############a 3 dimensional function I use for testing
detpep10curv <- function(xx){
  ##########################################################################
  #
  # DETTE & PEPELYSHEV (2010) CURVED FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUT:
  #
  # xx = c(x1, x2, x3)
  #
  #########################################################################
  
  x1 <- xx[1]
  x2 <- xx[2]
  x3 <- xx[3]
  
  term1 <- 4 * (x1 - 2 + 8*x2 - 8*x2^2)^2
  term2 <- (3 - 4*x2)^2
  term3 <- 16 * sqrt(x3+1) * (2*x3-1)^2
   y <- term1 + term2 + term3               ###no noise
 # y <- term1 + term2 + term3+rnorm(1,0,.5) ###noisy one
  return(y)
}

borehole <- function(xx)
{
  ##########################################################################
  #
  # BOREHOLE FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # OUTPUT AND INPUT:
  #
  # y  = water flow rate
  # xx = c(rw, r, Tu, Hu, Tl, Hl, L, Kw)
  #
  ##########################################################################
  
  rw <- xx[1]
  r  <- xx[2]
  Tu <- xx[3]
  Hu <- xx[4]
  Tl <- xx[5]
  Hl <- xx[6]
  L  <- xx[7]
  Kw <- xx[8]
  
  frac1 <- 2 * pi * Tu * (Hu-Hl)
  
  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1+frac2a+frac2b)
  
  y <- frac1 / frac2
  return(y)
}


fried <- function(xx)
{
  ##########################################################################
  #
  # FRIEDMAN FUNCTION
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUT:
  #
  # xx = c(x1, x2, x3, x4, x5)
  #
  ##########################################################################
  
  x1 <- xx[1]
  x2 <- xx[2]
  x3 <- xx[3]
  x4 <- xx[4]
  x5 <- xx[5]
  
  term1 <- 10 * sin(pi*x1*x2)
  term2 <- 20 * (x3-0.5)^2
  term3 <- 10*x4
  term4 <- 5*x5
  
  y <- term1 + term2 + term3 + term4
  return(y)
}
