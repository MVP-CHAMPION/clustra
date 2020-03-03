## This is the R translation of trajectories v2.sas code
## 
## This is multithreaded code that expects Unix OS with fork and 
## multithreaded OpenBLAS matrix support. Do not run on Mac or Win!
## 
## 
# library(data.table) # possibly for reading data
library(gratia)
library(parallel)
library(mgcv)
#library(proftools)
library(openblasctl)


## Function to fit thin plate spline (tps) to a group with gam from the mgcv
## package. Crossvalidation does not know about zero weights, resulting in 
## different smoothing parameters, so subset is used and a separate prediction
## function.
## TODO Add write generated data to a csv file
## TODO Add Ward's hierarchical clustering of the cluster means (evaluated on a
##      grid) for a converged high-K trajectories result.
##      Profile plots can help determine K.
## TODO Add uncertainty moderated outlier detection. Do we need to use credible
##      regions or is se.fit enough??
tps_g = function(g, dat, maxdf, nthreads) {
  tps = mgcv::bam(response ~ s(time, k = maxdf), data = dat, 
                  subset = dat$group == g, discrete = TRUE, nthreads = nthreads)
  if(class(tps)[1] == "try-error") print(attr(tps, "condition")$message)
  tps
}

## Function to predict for new data based on fitted tps of a group
pred_g = function(g, tps, data)
  predict(object = tps[[g]], newdata = data, type = "response", se.fit = TRUE)
## Loss functions for each id to a group tps center
mse_g = function(g, pred, dat) # mean squared error
  as.numeric(tapply((dat$response - pred[[g]]$fit)^2, dat$id, mean))
mxe_g = function(g, tps, dat) # maximum error
  as.numeric(tapply(abs(dat$response - pred[[g]]$fit), dat$id, max))
## Fn to compute points outside CI
## Mean distance in probability (a multi-Mahalanobis distance)
mpd_g = function(g, pred, dat) {
  ## TODO 
}

## Function to assign starting groups. 
#' @param data Data frame with response measurements, one per observation. Column
#' names are id, time, response, group. Note that id are already sequential
#' starting from 1. This affects expanding group numbers to ids.
#' @param ng Number of clusters (groups).
#' "bestrand" is implemented.
#' @param starts The number of random start values to consider.
#' @param start_nid The number of id's per cluster to sample for evaluating starts.
#' @param maxdf Maximum degrees of freedom for trajectory spline smooths.
#' @export
start_groups = function(data, ng, starts, start_nid, cores, maxdf) {
  cat("\nStarts : ")
  ## test data for diversity criterion
  test_data = data.frame(id = rep(0, ng*2*maxdf),
                time = rep(seq(min(data$time), max(data$time),
                               length.out = 2*maxdf), times = ng),
                response = rep(NA, ng*2*maxdf),
                group = rep(1:ng, each = 2*maxdf))

  ## Samples id's. Picks tps fit with best deviance after one iteration among
  ##   random starts. Choosing from samples increases diversity of 
  ##   fits (sum of distances between group fits).
  ## Then classifies all ids based on fit from best sample
  max_div = 0
  start_id = sort(sample(unique(data$id), ng*start_nid)) # sample to add diversity!
  dat_start = dat[data$id %in% start_id, ]
  dat_start$id = as.numeric(factor(dat_start$id)) # since id are not sequential (Am I losing correspondence to ids? No, because spline modesl do not know about ids. The classification process only needs who has same id, not what is the id.)
  for(i in 1:starts) {
    group = sample(ng, ng*start_nid, replace = TRUE) # random groups
    dat_start$group = group[dat_start$id] # expand group to all responses 
    
    f = trajectories(dat_start, ng, group, iter = 1, maxdf = maxdf,
                     plot = FALSE, cores = cores)
    if(any(lapply(f$tps, class) == "try-error")) next
#    deviance = sum(unlist(lapply(1:ng, function(g) deviance(f$tps[[g]]))))
   
    diversity = sum(dist(
      do.call(rbind, lapply(mclapply(1:ng, pred_g, tps = f$tps, data = test_data,
                              mc.cores = cores$e_mc), function(x) as.numeric(x$fit)))))
    if(diversity > max_div) {
      max_div = diversity
      best_tps_cov = f$tps
    }
    cat(round(diversity, 2), "")
  }
  cat("->", max_div, "")
  ## predict for all observations of all ids
  pred = mclapply(1:ng, pred_g, tps = best_tps_cov, data = data,
                  mc.cores = cores$e_mc)
  ## compute mse for each id
  loss = mclapply(1:ng, mse_g, pred = pred, dat = data, 
                  mc.cores = cores$e_mc) # !!! loss is not dat dimensioned!!!!
  ## classify id to group with min mse
  group = apply(do.call(cbind, loss), 1, which.min)
  group ## report group for each id
}


#' Cluster longitudinal trajectories of a response variable. Trajectory means
#' are splines fit to all ids in a cluster.
#' @param dat Data frame with response measurements, one per observation. Column
#' names are id, time, response, group. Note that id are already sequential
#' starting from 1. This affects expanding group numbers to ids.
#' @param ng Number of clusters (groups)
#' @param group Group numbers corresponding to sequential ids.
#' @param iter Maximum number of iterations. Note if iter <= 2, only deviance is
#' displayed for minimal output during start values exploration.
#' @param maxdf Maximum degrees of freedom for trajectory spline smooths.
#' @param plot Plot clustered data with superimposed trajectory spline smooths.
#' @param cores A list specifying multicore parallelism with 
#' components e_mc (expectation across ng), m_mc (maximization across ng), 
#' bam_nthreads (see bam documentation), blas (OpenBlas). Care should be taken
#' that cores are not oversubscribed.
#' @export
trajectories = function(dat, ng, group, iter = 15, maxdf = 50, plot = FALSE,
                        cores = list(e_mc = 1, m_mc = 1, bam_nthreads = 1, blas = 1)) {
  openblas_set_num_threads(cores$blas)
  if(max(dat$id) != length(group)) 
    cat("\ntrajectories: imput id's NOT sequential!\n")
  if(iter > 2) catlev = 1 else catlev = 0
  
  ## get number of unique id
  if(plot) require(ggplot2)
  n_id = length(unique(dat$id)) 
  
  ## EM algorithm to cluster ids into ng groups
  ## iterate fitting a thin plate spline center to each group (M-step)
  ##         regroup each id to nearest tps center (E-step)
  if(catlev) a = a_it = deltime(a)
  for(i in 1:iter) {
    if(catlev) cat("\n", i, "")
    ##
    ## M-step:
    ##   Estimate tps model parameters for each cluster from id's in that cluster
    tps = mclapply(1:ng, tps_g, dat = dat, maxdf = maxdf,
                   mc.cores = cores$m_mc, nthreads = cores$bam_nthreads)
    if(catlev) a = deltime(a, "M-step")

    if(any(lapply(tps, class) == "try-error")) {
      cat("\n")
      for(g in 1:ng) if(class(tps[[g]])[1] == "try-error") 
        print(paste0("Group ", g, " :", attr(tps[[g]], "condition")$message))
      cat("Random reshuffle for next iteration.\n")
      new_group = sample(ng, n_id, replace = TRUE)
      changes = sum(new_group != group)
      counts = tabulate(new_group)
    } else {
      ##
      ## E-step:
      ##   predict each id trajectory from each tps model
      pred = mclapply(1:ng, pred_g, tps = tps, dat = dat, mc.cores = cores$e_mc)
      ##   compute loss of each id to each to spline
      loss = mclapply(1:ng, mse_g, pred = pred, dat = dat, mc.cores = cores$e_mc)
      ## classify each id to mse-nearest tps model
      new_group = apply(do.call(cbind, loss), 1, which.min)
      if(catlev) a = deltime(a, " E-step")
      changes = sum(new_group != group)
      counts = tabulate(new_group)
      deviance = sum(unlist(lapply(1:ng, function(g) deviance(tps[[g]]))))
      if(catlev) 
         cat(" Changes:", changes, "Counts:", counts, "Deviance:", deviance)
    }    
    group = new_group
    dat$group = as.factor(group[dat$id]) # expand group to data frame
    if(changes == 0) break
    
    ## add outlier treatment as probability of belonging?
    ## use predict.gam(), probably for each group separately? weights??
  }
  if(catlev) a = deltime(a_it, " Done")
  
  if(plot) {
    plot_tps(dat)
    a = deltime(a, "\nDone plots")
  }
  list(deviance = deviance, group = group, dat_group = dat$group, tps = tps)
}

