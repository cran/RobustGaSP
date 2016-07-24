##########################################################################
## rgasp fit function
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
#library(nloptr)    ####I need to have this nloptr, or I can use C++ nlopt, this one is faster than optim()
#requireNamespace("nloptr")
                 
rgasp <- function(design, response,trend=matrix(1,length(response),1),nugget=0,
                  nugget.est=F,range.par=NA,prior_choice='ref_approx',a=0.2,
                  b=1/(length(response))^{1/dim(as.matrix(design))[2]}*(a+dim(as.matrix(design))[2]),
                  kernel_type='matern_5_2',
                  alpha=rep(1.9,dim(as.matrix(design))[2]),multiple_starts=T,lower_bound=T,max_eval=30,xtol_rel=1e-5){
  
  
  if (!is.logical(nugget.est) && length(nugget.est) != 1){
    stop("nugget.est should be boolean (either T or F) \n")
  }
  
  if(nugget!=0 & nugget.est==T){
    stop("one cannot fix and estimate the nugget at the same time \n")  
  }
  
  if(!is.na(range.par)){  
    if(length(range.par)!=dim(as.matrix(design))[2]){
      stop("range.par should either be fixed or estimated.")    
    }
    if(nugget.est){
      stop("We do not support fixing range parameters while estimating the nugget.")      
    }
  }
  
  if(!is.numeric(nugget)){
    stop("nugget should be a numerical value \n")  
  }
  
  model <- new("rgasp")
  model@call <- match.call()
  
  model@kernel_type <-kernel_type
  cat('Kernel type is ',model@kernel_type)
  
  model@alpha <-alpha
  
  model@input <- as.matrix(design)
  # print(model@input)
  model@output <- as.matrix(response)
  # print(model@output)
  model@p <- dim(model@input)[2]
  model@num_obs <- dim(model@input)[1]
  
  ## Checking the dimensions between the different inputs
  if (model@num_obs != dim(model@output)[1]){
    stop("The dimensions of the design matrix and the response do not match. \n")
  }
  #  stop("The dimension of the training points and simulator values are different. \n")
  
  #####I now build the gasp emulator here
  ##default setting constant mean basis. Of course We should let user to specify it
  #model@X = matrix(1,model@num_obs,1)   ####constant mean
  model@X=trend;               ###If the trend is specified, use it. If not, use a constant mean. 
  #######number of regressor
  model@q = dim(model@X)[2]; # NOTE THIS IS ALWAYS 1 SINCE YOU DEFINE IT THAT WAY ABOVE
  ####################correlation matrix
  
  model@nugget.est <- nugget.est
  
  model@R0 = as.list(1:model@p)
  for(i in 1:model@p){
    model@R0[[i]] = as.matrix(abs(outer(model@input[,i], model@input[,i], "-")))
  }
  
  ###########calculating lower bound for beta
  model@CL = rep(0,model@p)    ###CL is also used in the prior so I make it a model parameter
  
  for(i_cl in 1:model@p){
    model@CL[i_cl] = mean(model@R0[[i_cl]][which(model@R0[[i_cl]]>0)])
  }
  
  if(is.na(range.par)){
    ########this also depends on the kernel
    ##get some initial values of beta that we can start from a good set of range parameters to optimize
    COND_NUM_UB = 10^{16}  ###maximum condition number, this might be a little too large
    
    LB_all = optimize(search_LB_prob, interval=c(-5,12), maximum = FALSE, R0=model@R0,COND_NUM_UB= COND_NUM_UB,
                      p=model@p,kernel_type=model@kernel_type,alpha=model@alpha,nugget=nugget) ###find a lower bound for parameter beta
    
    LB_prob = exp(LB_all$minimum)/(exp(LB_all$minimum)+1)
    
    LB = NULL
    
    for( i_LB in 1:model@p){
      LB = c(LB, log(-log(LB_prob)/(max(model@R0[[i_LB]]))))    ###LB is lower bound for log beta, may consider to have it related to p
    }
    
    ##LB=LB-log(model@p)   #####

    if(lower_bound==T){
      if(model@nugget.est){
        model@LB=c(LB,-Inf)   
      }else{
        model@LB=LB
      }
    }else{
      if(model@nugget.est){
        model@LB=rep(-Inf,model@p+1)
      }else{
        model@LB=rep(-Inf,model@p)
      }
    }
    
    ############################the lower bound might be needed to discuss
    
    #beta_initial=(a+model@p)/(model@p*model@CL*b)/2  ###half the prior value
    beta_initial=50*exp(LB) #####one start
    
    model@beta_initial=beta_initial
    
    if(model@nugget.est){
      ini_value=c(log(beta_initial), log(0.0001)); ###with nugget
    }else{
      ini_value=log(beta_initial);               ###without the nugget
    }
    
    if(prior_choice=='ref_approx'){####this one can be with nugget or without the nugget
   #   if (requireNamespace("lbfgs", quietly = TRUE)) {
      tt_all <- nloptr::lbfgs(ini_value, neg_log_marginal_post_approx_ref, 
                      neg_log_marginal_post_approx_ref_deriv,nugget=nugget, nugget.est=model@nugget.est, 
                      R0=model@R0,X=model@X, output=model@output, CL=model@CL, a=a,b=b,
                      kernel_type=model@kernel_type,alpha=model@alpha,lower=model@LB,
                      nl.info = TRUE, control = list(maxeval=max_eval, xtol_rel=xtol_rel))
    #  }
    }else if(prior_choice=='ref_xi'|prior_choice=='ref_gamma'){####this needs to be revised
    #  if (requireNamespace("lbfgs", quietly = TRUE)) {
        tt_all <- nloptr::lbfgs(ini_value, neg_log_marginal_post_ref, 
                      nugget=nugget, nugget.est=nugget.est, R0=model@R0,
                      X=model@X, output=model@output, prior_choice=prior_choice, kernel_type=model@kernel_type,
                      alpha=model@alpha,lower=model@LB,nl.info = TRUE, control = list(maxeval=max_eval, xtol_rel=xtol_rel))
    #  }
    }

    ###########here I do multiple starts 
    if(multiple_starts){
      beta_initial=(a+model@p)/(model@p*model@CL*b)/2  ###half the prior value
      
      #  beta_initial=(a+model@p)/(model@p*model@CL*b)/2 #####another start
      
      if(nugget.est==T){
        ini_value=c(log(beta_initial), log(0.0002)); ###with nugget
        
      }else{
        ini_value=log(beta_initial);               ###without the nugget
      }
      
      if(prior_choice=='ref_approx'){####this one can be with nugget or without the nugget
      #  if (requireNamespace("lbfgs", quietly = TRUE)) {
          tt_all2 <- nloptr::lbfgs(ini_value, neg_log_marginal_post_approx_ref, 
                                  neg_log_marginal_post_approx_ref_deriv,nugget=nugget, nugget.est=nugget.est, R0=model@R0,
                                  X=model@X, output=model@output, CL=model@CL, a=a,b=b,  kernel_type=kernel_type,alpha=alpha,lower=model@LB,
                                  nl.info = TRUE, control = list(maxeval=max_eval, xtol_rel=xtol_rel))
      #  }
      }else if(prior_choice=='ref_xi'|prior_choice=='ref_gamma'){####this needs to be revised
       # if (requireNamespace("lbfgs", quietly = TRUE)) {
          tt_all2 <- nloptr::lbfgs(ini_value, neg_log_marginal_post_ref, 
                                  nugget=nugget, nugget.est=nugget.est, R0=model@R0,
                                  X=model@X, output=model@output, prior_choice=prior_choice, kernel_type=kernel_type,alpha=alpha,lower=model@LB,
                                  nl.info = TRUE, control = list(maxeval=max_eval, xtol_rel=xtol_rel))
      #  }
      }
      
      if(tt_all2$value<tt_all$value){
        model@beta_initial=beta_initial
        tt_all=tt_all2
      }
    }
    
    ################################
    ##############estimated beta
    # print out to see whether or not it has converged
   # cat('Convergence: ',tt_all,'\n')  ###cat doesn't work here
    
   
    model@log_post=-tt_all$value
    
    model@nugget=nugget;
    if(nugget.est){
      model@beta = exp(tt_all$par)[1:model@p];
      model@nugget=exp(tt_all$par)[model@p+1];
    }else{
      model@beta = exp(tt_all$par);
      #  model@nugget=0;
    }
  }else{###this is the case where the range parameters and the nugget are all fixed
    model@LB=NULL
    model@beta=range.par
    model@nugget=nugget
  }
  
  list_return=construct_rgasp(model@beta, model@nugget, model@R0, model@X, 
                              model@output,model@kernel_type,model@alpha); 
  model@L=list_return[[1]];
  model@LX=list_return[[2]];
  model@theta_hat=list_return[[3]];
  model@sigma2_hat=list_return[[4]];
  
  return(model)
}

show.rgasp <- function(object) {	
  cat("\n")
  cat("Call:\n")
  print(object@call)
  
  # cat('The dimension names of the design are: ',colnames(object@input),'\n')
  #  cat('The dimension names of the response is: ',colnames(object@output),'\n')
  cat('The dimension  of the design is: ',dim(object@input),'\n')
  cat('The number of the output is: ',dim(object@output),'\n')
  cat('The slots of this object are: ',slotNames(object),'\n')
  
}