## This is the R translation of trajectories v2.sas code
library(data.table)
library(tidyverse)

dir = "~/SciApplications/VA-MVP/03-Longitudinal/Trajectories"
set.seed(88)
setwd(dir)

source("generate.R")
id_df = gen_bp_data(n_id = 1000, m_obs = 25)

trajectories = function(id_df, ngroups = 5, iter = 10, maxdf = 50) {
  ## start with random group assignments
  df = id_df %>% group_by(id) %>% 
    mutate(group = rep(sample(ngroups, 1), length(id))) %>% 
    ungroup()
  
  ## spline fitting function
  ## gam from the mgcv package?? Includes crossvalidation.
  splineit = function()
  
  ## EM algorithm to cluster ids into ngroups
  ## iterate fitting a spline center to each group (M-step)
  ##         assigning each individual to nearest spline center (E-step)
  for(i in 1:iter) {
    ## M-step
    
  }
}


## Comments from DG:
## 
## doing same thing with Framingham data
##  bp data
##  id links records
##    1 record per bloodpressure
##    
##  with va one or more bp measurements
##  
##  Purpose is to end up with subjects are clustered with similarity of bp trajectories
##  
##  2 pieces: macro trajectories has everyt
##  
##  dsn = dataset name
##  id = 
##  time variable
##  risk variable
##  number of groups  ~ 5
##  max iterations ~ 50
##  maxdf - max percentage allowed to move
##  
##  assign people randomly to clusters - random starts
##  
##  do - number of iterations
##  
##  fitting a spline  gampl - generalized additive model spline
##    every subject is in each model but only ones in group have nonzero weight (1-0)
##    
##  high performance = 8 cores
##  
##  do this for each group
##  take output data for each group
##  
##  one obs per bp measurement - 5 columns of predicted values
##           take (predicted - actual)^2
##           take mean
##           
##  Also have a trim procedure to take out outliers
##  
##  Rest is just plotting
##  
##  Convert to R 
##  1. do it for multivariate
##  
##  do 100 clusters 
##  add george to KDI space - 
##  