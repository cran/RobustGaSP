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

   if( dim(testing_trend)[2]!=dim(object@X)[2]){
     stop("The dimensions of the design matrix and testing inputs matrix do not match. \n")
  
   }
#   ########################set up the mean and the correlation matrix
#   
#   if(!Ccall){
#   num_testing_input <- dim(testing_input)[1]
#   
#     
#   ####should let users specify the trend
#   X_testing = matrix(1,num_testing_input,1) 
#     ###################caculate correlation between the input and testing input
#   
#   r=matrix(0,num_testing_input,object@num_obs)######
#   r0=matrix(0,object@num_obs,object@dim_inputs);
#   
#   for(qi in  1: num_testing_input ){
#     
#     for(i in 1:(object@dim_inputs) ){     
#      # input_star=testing_input[qi,i];
#       r0[,i]=Matern_2_5_funct(abs(testing_input[qi,i]-object@input[,i]),object@beta[i])
#     }
#     
#     r[qi,]=r0[,1]
#     if(object@dim_inputs>1){
#       for(i in 2:object@dim_inputs ){
#         r[qi,]=r[qi,]*r0[,i]
#       }
#     }
#   }
#   
#   r_R_inv=t(backsolve(t(object@L),forwardsolve(object@L, t(r)) ))
#   ###############c_star_star for the variance in my Roubst GaSP paper
#   c_star_star=rep(0,num_testing_input)
#   
#   for(i in 1:num_testing_input){
#     X_testing_X_R.inv_r_i=(X_testing[i,]-r[i,]%*%object@R.inv_X );
#     c_star_star[i]=object@R[1,1]-r_R_inv[i,]%*%r[i,]+(X_testing_X_R.inv_r_i)%*%object@Xt_R.inv_X.inv%*%t(X_testing_X_R.inv_r_i)
#     
#   }
#   
#   
#   #########compute the predictive mean
#   ####same as MU_testing1  
#   MU_testing=  X_testing%*%object@theta_hat+(r_R_inv)%*%(object@output-object@X%*%object@theta_hat)        ####same as MU_testing1
#   
#   #########compute the lower and upper bound
#   
#  # var=matrix(c_star_star)%*%matrix(object@S2)/(object@num_obs-object@q) ###this hasn't take into account the variance of t
#   
#   var=matrix(c_star_star)%*%matrix(object@sigma2_hat) ###this hasn't take into account the variance of t
#  
#  
#   LB_testing=MU_testing+sqrt(var)*qt(0.025,df=(object@num_obs-object@q))
#   
#   UB_testing=MU_testing+sqrt(var)*qt(0.975,df=(object@num_obs-object@q))
#   
#   
#   
# 
#   output.list <- list()
#   output.list$mean=MU_testing   #####can we all use @ or S? It will be more user friendly in that way 
#   output.list$lower95=LB_testing
#   output.list$upper95=UB_testing
#   output.list$sd=var*(object@num_obs-object@q)/(object@num_obs-object@q-2)  ##the second part is t variance
 # }else{
    ####should let users specify the trend
    num_testing_input <- dim(testing_input)[1]
    #X_testing = matrix(1,num_testing_input,1) ###testing trend
      
    qt_025=qt(0.025,df=(object@num_obs-object@q))
    qt_975=qt(0.975,df=(object@num_obs-object@q))
      
      
    testing_input=as.matrix(testing_input)

    r0 = as.list(1:object@p)
    for(i in 1:object@p){
      r0[[i]] = as.matrix(abs(outer(testing_input[,i], object@input[,i], "-")))
    }
      

    pred_list=pred_rgasp(object@beta,object@nugget,object@input,object@X,object@output,
                         testing_input,testing_trend,object@L,object@LX,object@theta_hat,
                         object@sigma2_hat,qt_025,qt_975,r0,object@kernel_type,object@alpha)
    
    output.list <- list()
    output.list$mean=pred_list[[1]]   #####can we all use @ or S? It will be more user friendly in that way 
    output.list$lower95=pred_list[[2]]
    output.list$upper95=pred_list[[3]]
    output.list$sd=pred_list[[4]] 
#  }
 
#     
#   output.list$MSE_emulator <- MSE_emulator
#   output.list$prop_emulator <- prop_emulator
#   output.list$length_emulator <- length_emulator
#   output.list$cond_num_R <- cond_num_R
  
  return (output.list)
}
