## This is the R translation of trajectories v2.sas code
library(data.table)
library(tidyverse)
library(mgcv)
library(parallel)

dir = "~/Git/mvp-champion/trajectories/"
set.seed(89)
setwd(dir)

source("R/generate.R")
n_id = 1000
df = gen_bp_data(n_id = n_id, m_obs = 25, plots = 10)

trajectories = function(df, ngroups = 3, iter = 15, maxdf = 50) {
  ## FIXME n_id as parameter
  ## start with random group assignments
  group = sample(ngroups, n_id, replace = TRUE)
  df_id = df$id
  df_group = group[df_id]
  
  ## Fn to fit thin plate spline (tps) to a group with gam from mgcv package
  spline_k = function(k, df_group) 
    gam(sbp ~ s(age, k = maxdf), data = df, weights = df_group == k)
  ## Fn to compute mse for each id to a tps center
  mse_k = function(x, fit) 
    tapply((df$sbp - fit[[x]]$fitted.values)^2, INDEX = df$id, FUN = mean)

  ## EM algorithm to cluster ids into ngroups
  ## iterate fitting a spline center to each group (M-step)
  ##         regroup each id to nearest spline center (E-step)
  for(i in 1:iter) {
    ## M-step:
    ##   fit tp spline centers for each group
    tps_mean = mclapply(1:ngroups, FUN = spline_k, df_group = df_group, mc.cores = 2)
    
    ## E-step:
    ##   compute mse to each tp spline
    mse = mclapply(1:ngroups, FUN = mse_k, fit = tps_mean, mc.cores = 2)
    ## get nearest tp spline
    new_group = apply(do.call(cbind, mse), 1, which.min)
    
    ## done?
    changes = sum(new_group != group)
    cat("iteration", i, "changes", changes, "\n")
    group = new_group
    df_group = group[df_id]
    if(changes == 0) break
    
    ## add outlier treatment
  }
  group
}
system.time((df = trajectories(id_df)))
