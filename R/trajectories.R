## This is the R translation of trajectories v2.sas code
library(data.table)
library(tidyverse)
library(mgcv)

dir = "~/SciApplications/VA-MVP/03-Longitudinal/Trajectories"
set.seed(88)
setwd(dir)

source("generate.R")
id_df = gen_bp_data(n_id = 1000, m_obs = 25)

trajectories = function(id_df, ngroups = 5, iter = 1, maxdf = 50) {
  ## start with random group assignments
  df = id_df %>% group_by(id) %>% 
    mutate(group = rep(sample(ngroups, 1), length(id))) %>% 
    ungroup()
  
  ## spline fitting function
  ## gam from the mgcv package?? Includes crossvalidation.
  ##   subset parameter or by variable faster?
  spline_k = function(k) gam(sbp ~ s(age, k = maxdf),
                             data = filter(df, group == k),
                             xeval = filter(df, group != k))
  
  ## EM algorithm to cluster ids into ngroups
  ## iterate fitting a spline center to each group (M-step)
  ##         assigning each individual to nearest spline center (E-step)
  for(i in 1:iter) {
    ## M-step: fit tp
    f_mean = lapply(1:ngroups, spline_k)
    
    ## E-step
  }
  f_mean
}
f_mean = trajectories(id_df)
