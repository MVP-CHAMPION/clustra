## This is the R translation of trajectories v2.sas code
library(data.table)
library(tidyverse)
library(mgcViz)
library(parallel)
library(pbdIO)
a0 = a = deltime()

#' Cluster longitudinal trajectories of a response variable. Trajectory means
#' are splines fit to all ids in a cluster.
#' @param dat Data frame with response measurements, one per observation. Column
#' names are id, time, response.
#' @param ngroups 
#' @param iter Maximum number of iterations.
#' @param maxdf Maximum degrees of freedom for trajectory spline centers.
#' @export
trajectories = function(dat, ngroups, iter = 20, maxdf = 50) {
  ## get number of unique id
  n_id = length(unique(dat$id)) 
  
  ## start with random group assignments
  group = sample(ngroups, n_id, replace = TRUE)
  dat_id = dat$id
  dat_group = group[dat_id] # expand group assignment to all responses
  
  ## Function to fit thin plate spline (tps) to a group with gam from mgcv package
  spline_k = function(k, dat_group, dat, maxdf) 
    mgcv::gam(response ~ s(time, k = maxdf), data = dat, weights = dat_group == k)
  ## Loss function for each id to a tps center
  mse_k = function(x, tps_mean, dat_id) # mean squared error
    tapply((dat$response - fitted(tps_mean[[x]]))^2, INDEX = dat_id, FUN = mean)
  mxe_k = function(x, tps_mean, dat_id) # maximum error
    tapply(abs(dat$response - fitted(tps_mean[[x]])), INDEX = dat_id, FUN = max)
  ## Fn to compute points outside CI

  ## EM algorithm to cluster ids into ngroups
  ## iterate fitting a thin plate spline center to each group (M-step)
  ##         regroup each id to nearest spline center (E-step)
  a = deltime(a, "Starting iteration")
  for(i in 1:iter) {
    ## M-step:
    ##   fit tp spline centers for each group separately (via weights of gam)
    tps_mean = mclapply(1:ngroups, FUN = spline_k, dat_group = dat_group, dat = dat,
                        maxdf = maxdf, mc.cores = 4)
    a = deltime(a, "M-step")
    ## E-step:
    ##   compute MSE of each id to each tp spline
    mse = sapply(1:ngroups, FUN = mse_k, tps_mean = tps_mean, dat_id = dat_id)
    ## get mse-nearest tp spline for each id
    new_group = apply(mse, 1, which.min)
    a = deltime(a, "E-step")
    
    ## Done?
    changes = sum(new_group != group)
    cat("iteration", i, "changes", changes, "\n")
    group = new_group
    dat_group = group[dat_id]
    if(changes == 0) break
    
    ## add outlier treatment as probability of belonging?
      
  }
  list(group = group, tps_mean = tps_mean)
}

source("R/generate.R")
setwd("~/Git/mvp-champion/trajectories/")
set.seed(8837)

dat = gen_long_data(n_id = 1000, m_obs = 25, plots = 20)
a = deltime(a, "Data generated")

## Rprof() # uncomment to profile time
gam_ctl = gam.control(nthreads = 1)
f = trajectories(dat, 3)
## Rprof(NULL)
## summaryRprof()

g1 = getViz(f$tps_mean[[1]])
g2 = getViz(f$tps_mean[[2]])
g3 = getViz(f$tps_mean[[3]])
p1 = plot( sm(g1, 1) ) + l_fitLine(colour = "red") +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()
p2 = plot( sm(g2, 1) ) + l_fitLine(colour = "red") +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()
p3 = plot( sm(g3, 1) ) + l_fitLine(colour = "red") +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()
print(p1)
print(p2)
print(p3)

a = deltime(a0, "Total time")

