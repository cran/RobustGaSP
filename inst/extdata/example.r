##########################################################################
## example.R
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

#############
# Simulated example
#############

library("RobustGaSP")

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
# m1<- rgasp(design = input, response = output, lower_bound=T)

# the following use a single start on optimization
# m1<- rgasp(design = input, response = output, lower_bound=F, multiple_starts=F)

# number of points to be predicted 
num_testing_input <- 5000    
# generate points to be predicted
testing_input <- matrix(runif(num_testing_input*dim_inputs),num_testing_input,dim_inputs)
# Perform prediction
m1.predict<-predict(m1, testing_input)
# Predictive mean
m1.predict$mean  

# The following tests how good the prediction is 
testing_output <- matrix(0,num_testing_input,1)
for(i in 1:num_testing_input){
  testing_output[i]<-detpep10curv(testing_input[i,])
}

# compute the MSE, average coverage and average length
# out of sample MSE
MSE_emulator <- sum((m1.predict$mean-testing_output)^2)/(num_testing_input)  

# proportion covered by 95\% posterior predictive credible interval
prop_emulator <- length(which((m1.predict$lower95<=testing_output)
                              &(m1.predict$upper95>=testing_output)))/num_testing_input

# average length of  posterior predictive credible interval
length_emulator <- sum(m1.predict$upper95-m1.predict$lower95)/num_testing_input

# output of prediction
MSE_emulator
prop_emulator
length_emulator  
# normalized RMSE
sqrt(MSE_emulator/mean((testing_output-mean(output))^2 ))
