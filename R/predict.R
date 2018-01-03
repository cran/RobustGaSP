##########################################################################
## prediction function
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

predict.rgasp <- function (object, testing_input, testing_trend= matrix(1,dim(testing_input)[1],1),...){
  if(object@zero_mean=="Yes"){
    testing_trend=rep(0,dim(testing_input)[1]);
  }else{
    if( dim(testing_trend)[2]!=dim(object@X)[2]){
      stop("The dimensions of the design trend matrix and testing trend matrix do not match. \n")
    }
  }
  
  if( dim(testing_input)[2]!=dim(object@input)[2]){
    stop("The dimensions of the design matrix and testing inputs matrix do not match. \n")
  }
  ####let users specify the trend
    num_testing_input <- dim(testing_input)[1]
    #X_testing = matrix(1,num_testing_input,1) ###testing trend
      
    qt_025=qt(0.025,df=(object@num_obs-object@q))
    qt_975=qt(0.975,df=(object@num_obs-object@q))
      
      
    testing_input=as.matrix(testing_input)

    r0 = as.list(1:object@p)
    for(i in 1:object@p){
      r0[[i]] = as.matrix(abs(outer(testing_input[,i], object@input[,i], "-")))
    }
      

    pred_list=pred_rgasp(object@beta_hat,object@nugget,object@input,object@X,object@zero_mean,object@output,
                         testing_input,testing_trend,object@L,object@LX,object@theta_hat,
                         object@sigma2_hat,qt_025,qt_975,r0,object@kernel_type,object@alpha)
    
    output.list <- list()
    output.list$mean=pred_list[[1]]   #####can we all use @ or S? It will be more user friendly in that way 
    output.list$lower95=pred_list[[2]]
    output.list$upper95=pred_list[[3]]
    output.list$sd=sqrt(pred_list[[4]]) 

  return (output.list)
}
