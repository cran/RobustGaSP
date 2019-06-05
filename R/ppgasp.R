##########################################################################
## rgasp fit function
## 
## Robust GaSP Package
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, April 2013.
##
## Copyright (C) 2015-present Mengyang Gu, Jesus Palomo, James O. Berger
##							  
##    
##########################################################################

##function for parallel partial Gaussian stochastic process
ppgasp <- function(design, response,trend=matrix(1,dim(response)[1],1),zero.mean="No",nugget=0,
                  nugget.est=F,range.par=NA,prior_choice='ref_approx',a=0.2,
                  b=1/(length(response))^{1/dim(as.matrix(design))[2]}*(a+dim(as.matrix(design))[2]),
                  kernel_type='matern_5_2',
                  alpha=rep(1.9,dim(as.matrix(design))[2]),lower_bound=T,max_eval=max(30,20+5*dim(design)[2]),
                  initial_values=NA,num_initial_values=2){
  
  ##response is n x k matrix, where n is the number of obs and k is the number of loc
  if (zero.mean=="Yes"){
    trend=matrix(0,dim(response)[1],1)
  }
  
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
  
  #model <- new("rgasp")
  model <- new("ppgasp")
  
  model@call <- match.call()
  
  
  #cat('Kernel type is ',model@kernel_type,'\n')
  
  #cat('Multiple starting points: ',multiple_starts, '\n')
  
  model@alpha <-alpha
  
  #####Only numeric inputs are allowed
  design=as.matrix(design)
  model@input <- matrix(as.numeric(design), dim(design)[1],dim(design)[2])
  
  # print(model@input)
  model@output <- as.matrix(response)
  
  if(dim(model@output)[1]!=dim(model@input)[1]){
    stop("the number of rows in the input should be the same as the number of rows 
         in the output \n")  
    
  }
  # print(model@output)
  model@p <- dim(model@input)[2]
  model@num_obs <- dim(model@output)[1]
  model@k <- dim(model@output)[2]
  
  ## Checking the dimensions between the different inputs
  if (model@num_obs != dim(model@output)[1]){
    stop("The dimensions of the design matrix and the response do not match. \n")
  }
  #  stop("The dimension of the training points and simulator values are different. \n")
  

  ##kernel type
  if(length(kernel_type)==1){
    kernel_type=rep(kernel_type,model@p)
  }else if(length(kernel_type)!=model@p){
    stop("Please specify the correct number of kernels. \n")
  }
  model@kernel_type <-kernel_type
  
  ##change kernel type to integer to pass to C++ code
  ##1 is power exponent, 2 is matern with roughness 3/2, and 3 is matern with roughenss parameter 5/2
  kernel_type_num=rep(0,  model@p)
  for(i_p in 1:  model@p){
    if(model@kernel_type[i_p]=="matern_5_2"){
       kernel_type_num[i_p]=as.integer(3)
    }else if (model@kernel_type[i_p]=="matern_3_2"){
      kernel_type_num[i_p]=as.integer(2)
    }else if (model@kernel_type[i_p]=="pow_exp"){
      kernel_type_num[i_p]=as.integer(1)
    }
  }
    


  #####I now build the gasp emulator here
  ##default setting constant mean basis. Of course We should let user to specify it
  #model@X = matrix(1,model@num_obs,1)   ####constant mean
  model@X=trend;               ###If the trend is specified, use it. If not, use a constant mean. 
  model@zero_mean=zero.mean
  #######number of regressor
  if(model@zero_mean=="Yes"){
    model@q=as.integer(0)
  }else{
    model@q = dim(model@X)[2]; # NOTE THIS IS ALWAYS 1 SINCE YOU DEFINE IT THAT WAY ABOVE
  }
  ####################correlation matrix
  
  model@nugget.est <- nugget.est
  
  model@R0 = as.list(1:model@p)
  for(i in 1:model@p){
    model@R0[[i]] = as.matrix(abs(outer(model@input[,i], model@input[,i], "-")))
  }
  
  ###########calculating lower bound for beta
  model@CL = rep(0,model@p)    ###CL is also used in the prior so I make it a model parameter
  
  for(i_cl in 1:model@p){
    #  model@CL[i_cl] = mean(model@R0[[i_cl]][which(model@R0[[i_cl]]>0)])
    
    model@CL[i_cl] = (max(model@input[,i_cl])-min(model@input[,i_cl]))/model@num_obs^{1/model@p}
  }
  
  
  if(is.na(range.par)){
    ########this also depends on the kernel
    ##get some initial values of beta that we can start from a good set of range parameters to optimize
    
   
      COND_NUM_UB = 10^{16}  ###maximum condition number, this might be a little too large
      
      LB_all = optimize(search_LB_prob, interval=c(-5,12), maximum = FALSE, R0=model@R0,COND_NUM_UB= COND_NUM_UB,
                        p=model@p,kernel_type=kernel_type_num,alpha=model@alpha,nugget=nugget) ###find a lower bound for parameter beta
      
      LB_prob = exp(LB_all$minimum)/(exp(LB_all$minimum)+1)
      
      LB = NULL
      
      for( i_LB in 1:model@p){
        LB = c(LB, log(-log(LB_prob)/(max(model@R0[[i_LB]]))))    ###LB is lower bound for log beta, may consider to have it related to p
      }
   
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
    
    cat('The upper bounds of the range parameters are',1/exp(model@LB),'\n')
    
    ############################the lower bound might be needed to discuss
    
    #beta_initial=(a+model@p)/(model@p*model@CL*b)/2  ###half the prior value
    
    if(is.na(initial_values)[1]==T){
      beta_initial=matrix(0,num_initial_values,model@p)
      eta_initial=rep(0,num_initial_values)
      beta_initial[1,]=50*exp(LB) #####one start
      eta_initial[1]=0.0001
      if(num_initial_values>1){
        beta_initial[2,]=(a+model@p)/(model@p*model@CL*b)/2  ###half the prior value
        eta_initial[2]=0.0002

      }
      if(num_initial_values>2){
        for(i_ini in 3:num_initial_values){
          set.seed(i_ini)
          
          beta_initial[i_ini,]=10^3*runif(model@p)/model@CL
          eta_initial[i_ini]=10^(-3)*runif(1)
        }
      }
      initial_values=cbind(log(beta_initial),log(eta_initial))
    }
      
    model@log_post=-Inf;
    for(i_ini in 1:num_initial_values){
        if(model@nugget.est){
          ini_value=initial_values[i_ini,]
        }else{
          ini_value=initial_values[i_ini,1:model@p]
          ###without the nugget
        }
        cat('The initial values of range parameters are', 1/exp(ini_value[1:model@p]),'\n')
        cat('Start of the optimization ', i_ini,' : \n')
        
        if(prior_choice=='ref_approx'){####this one can be with nugget or without the nugget
          #  if (requireNamespace("lbfgs", quietly = TRUE)) {
          tt_all <- try(nloptr::lbfgs(ini_value, neg_log_marginal_post_approx_ref_ppgasp, 
                                      neg_log_marginal_post_approx_ref_deriv_ppgasp,nugget=nugget, nugget.est=model@nugget.est, 
                                      R0=model@R0,X=model@X, zero_mean=model@zero_mean,output=model@output, CL=model@CL, a=a,b=b,
                                      kernel_type=kernel_type_num,alpha=model@alpha,lower=model@LB,
                                      nl.info = FALSE, control = list(maxeval=max_eval)),TRUE)
          #   }
        }else if(prior_choice=='ref_xi'|prior_choice=='ref_gamma'){####this needs to be revised
          #  if (requireNamespace("lbfgs", quietly = TRUE)) {
          tt_all <- try(nloptr::lbfgs(ini_value, neg_log_marginal_post_ref_ppgasp, 
                                      nugget=nugget, nugget.est=nugget.est, R0=model@R0,
                                      X=model@X, zero_mean=model@zero_mean,output=model@output, prior_choice=prior_choice, kernel_type=kernel_type_num,
                                      alpha=model@alpha,lower=model@LB,nl.info = FALSE, control = list(maxeval=max_eval)),TRUE)
          # }
        }
        #if(class(tt_all)=="try-error"){
          #sink()
        #  stop(tt_all)
        #}

        if(class(tt_all)!="try-error"){
          if(model@nugget.est==F){
            nugget_par=nugget
          }else{
            nugget_par=exp(tt_all$par)[model@p+1]
          }
          
          if (tt_all$convergence>0){convergence=T}else{convergence=F}
          cat('The number of interation is ', tt_all$iter,'\n',
              'The value of the posterior is ', -tt_all$value,'\n',
              'Optimized range parameters are', 1/exp(tt_all$par)[1:model@p],'\n',
              'Optimized nugget parameter is', nugget_par,'\n',
              'Convergence: ', convergence,'\n' )
        

          if( (-tt_all$value)>model@log_post){
            log_lik=-tt_all$value
            model@log_post=-tt_all$value
            model@nugget=nugget;
            if(nugget.est){
              model@beta_hat = exp(tt_all$par)[1:model@p];
              model@nugget=exp(tt_all$par)[model@p+1];
            }else{
              model@beta_hat = exp(tt_all$par);
              #  model@nugget=0;
            }
          }
        }
    }
  }else{###this is the case where the range parameters and the nugget are all fixed
    model@LB=rep(-Inf,model@p)
    model@beta_hat=1/range.par
    model@nugget=nugget
  }
  
  list_return=construct_ppgasp(model@beta_hat, model@nugget, model@R0, model@X, zero_mean=model@zero_mean,
                              model@output,kernel_type_num,model@alpha); 
  model@L=list_return[[1]];
  model@LX=list_return[[2]];
  model@theta_hat=list_return[[3]];
  model@sigma2_hat=list_return[[4]];
  
  return(model)
}

show.ppgasp <- function(object) {	
  cat("\n")
  cat("Call:\n")
  print(object@call)
  

  cat('The dimension  of the design is: ',dim(object@input),'\n')
  cat('The dimension  of the output is: ',dim(object@output),'\n')


  cat('Range parameters: ', 1/object@beta_hat,'\n')
  cat('Nugget parameter: ', object@nugget,'\n')
  
  
  # cat('The slots of this object are: ',slotNames(object),'\n')
  
}