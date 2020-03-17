## This is the R translation of trajectories v2.sas code
##
## This is multithreaded code that expects a Unix OS with fork (and
## multithreaded OpenBLAS matrix support - soon). Do not run on Mac or Win!
##
##

#' Function to set and maintain clustra parameters and working directory.
#' 
#' @details
#' Objective is to simplify clustra function invocation and manage sets of
#' clustra runs using a \code{jsonlite} file. Sets working directory to playdir.
#' Parameters PL set in the following order: if PL is supplied, it is written to
#' a json file and returned; if PL is NULL, attempts to find parname file in
#' playdir and use that; if PL is NULL and parname file is not present, PL is
#' set to defaults.
#' 
#' @param parname
#' Character string filename.
#' @param playdir 
#' Character string directory where filename is expected and written. 
#' @param PL 
#' A list of lists parameter data structure.
#' @export
clustra_par = function(parname = "clustra.par", playdir = "~/clustra_play",
                       PL = NULL) {
  setwd(playdir)
  if(!is.null(PL)) { # have new PL so write it out and use it
    jsonlite::write_json(PL, parname, pretty = TRUE)
  } else { # use PL from file or set default and write it
    if(file.exists(parname)) {
      PL = jsonlite::read_json(parname, simplifyVector = TRUE)
    } else {
      PL = list(
        gen_par = list(seed = 90, 
                       n_id = 10000, 
                       lambda_obs = 25,
                       first = c(-365, -10),
                       last = c(365*3, 365*10), 
                       plots = FALSE),
        cores = list(e_mc = 4, m_mc = 4, bam_nthreads = 1, blas = 1),
        traj_par = list(seed = 79,
                        maxdf = 30,
                        iter = 8, 
                        starts = list(ns = 5, nid = 10),
                        k_vec = c(2, 3, 4),
                        replicates = 4)
      )
      jsonlite::write_json(PL, parname, pretty = TRUE)
    }
  } 
  PL
}

library(gratia)
library(mgcv)
## TODO add controls to use one core if OpenBLAS is not available

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
tps_g = function(g, data, maxdf, nthreads) {
  tps = mgcv::bam(response ~ s(time, k = maxdf), data = data,
                  subset = data$group == g, discrete = TRUE, nthreads = nthreads)
  if(class(tps)[1] == "try-error") print(attr(tps, "condition")$message)
  tps
}

## Function to predict for new data based on fitted tps of a group
pred_g = function(g, tps, data)
  predict(object = tps[[g]], newdata = data, type = "response", se.fit = TRUE)
## Loss functions for each id to a group tps center
mse_g = function(g, pred, data) # mean squared error
  as.numeric(tapply((data$response - pred[[g]]$fit)^2, data$id, mean))
mxe_g = function(g, tps, data) # maximum error
  as.numeric(tapply(abs(data$response - pred[[g]]$fit), data$id, max))
## Fn to compute points outside CI
## Mean distance in probability (a multi-Mahalanobis distance)
mpd_g = function(g, pred, data) {
  ## TODO
}

#' Function to assign starting groups.
#' 
#' @param data Data frame with response measurements, one per observation.
#' Column names are id, time, response, group. Note that id are already
#' sequential starting from 1. This affects expanding group numbers to ids.
#' @param k Number of clusters (groups).
#' "bestrand" is implemented.
#' @param starts The number of random start values to consider.
#' @param snid The number of id's per cluster to sample in evaluating starts.
#' @param maxdf Maximum degrees of freedom for trajectory spline smooths.
#' @export
start_groups = function(data, k, nstart, nid, cores, maxdf, verbose = FALSE) {
  if(verbose) cat("\nStarts : ")
  ## test data for diversity criterion
  test_data = data.frame(id = rep(0, k*2*maxdf),
                time = rep(seq(min(data$time), max(data$time),
                               length.out = 2*maxdf), times = k),
                response = rep(NA, k*2*maxdf),
                group = rep(1:k, each = 2*maxdf))

  ## Samples k*nid id's. Picks tps fit with best deviance after one iteration among
  ##   random starts. Choosing from samples increases diversity of
  ##   fits (sum of distances between group fits).
  ## Then classifies all ids based on fit from best sample
  max_div = 0
  start_id = sort(sample(unique(data$id), k*nid)) # for speed and diversity!
  data_start = data[data$id %in% start_id, ]
  data_start$id = as.numeric(factor(data_start$id)) # since id are not sequential (Am I losing correspondence to ids? No, because spline modesl do not know about ids. The classification process only needs who has same id, not what is the id.)
  for(i in 1:nstart) {
    group = sample(k, k*nid, replace = TRUE) # random groups
    data_start$group = group[data_start$id] # expand group to all responses

    f = trajectories(data_start, k, group, iter = 1, maxdf = maxdf,
                     plot = FALSE, cores = cores)
    if(any(lapply(f$tps, class) == "try-error")) next

    diversity = sum(dist(
        do.call(rbind, lapply(mclapply(1:k, pred_g, tps = f$tps,
                                       data = test_data, mc.cores = cores$e_mc),
                              function(x) as.numeric(x$fit)))
    ))
    if(diversity > max_div) {
      max_div = diversity
      best_tps_cov = f$tps
    }
    if(verbose) cat(round(diversity, 2), "")
  }
  if(verbose) cat("->", max_div, "")
  ## predict for all observations of all ids
  pred = mclapply(1:k, pred_g, tps = best_tps_cov, data = data,
                  mc.cores = cores$e_mc)
  ## compute mse for each id
  loss = mclapply(1:k, mse_g, pred = pred, data = data,
                  mc.cores = cores$e_mc) # !!! loss is not data dimensioned!!!!
  ## classify id to group with min mse
  group = apply(do.call(cbind, loss), 1, which.min)
  group ## report group for each id
}


#' Cluster longitudinal trajectories of a response variable.
#' 
#' Trajectory means are splines fit to all ids in a cluster. Typically, this
#' function is called by \code{clustra()}.
#' 
#' @param datax Data frame with response measurements, one per observation. Column
#' names are \code{id}, \code{time}, \code{response}, \code{group}. Note that
#'  id are already sequential.
#' starting from 1. This affects expanding group numbers to ids.
#' @param k Number of clusters (groups)
#' @param group Group numbers corresponding to sequential ids.
#' @param iter Maximum number of iterations. Note if iter <= 2, only deviance is
#' displayed for minimal output during start values exploration.
#' @param maxdf Maximum degrees of freedom for trajectory spline smooths.
#' @param plot Plot clustered data with superimposed trajectory spline smooths.
#' @param cores A list specifying multicore parallelism with
#' components e_mc (expectation across k), m_mc (maximization across k),
#' bam_nthreads (see bam documentation), blas (OpenBlas). Care should be taken
#' that cores are not oversubscribed.
#' @export
trajectories = function(data, k, group, iter = 15, maxdf = 50, plot = FALSE,
                        cores = list(e_mc = 1, m_mc = 1, bam_nthreads = 1, blas = 1),
                        verbose = FALSE) {
  if(verbose) a = a_0 = deltime(a)
  library(parallel)
  openblasctl::openblas_set_num_threads(cores$blas)
  if(max(data$id) != length(group))
    cat("\ntrajectories: imput id's NOT sequential!\n")

  ## get number of unique id
  if(plot) require(ggplot2)
  n_id = length(unique(data$id))

  ## EM algorithm to cluster ids into k groups
  ## iterate fitting a thin plate spline center to each group (M-step)
  ##         regroup each id to nearest tps center (E-step)
  try_errors = 0
  for(i in 1:iter) {
    if(verbose) cat("\n", i, "")
    ##
    ## M-step:
    ##   Estimate tps model parameters for each cluster from id's in that cluster
    tps = mclapply(1:k, tps_g, data = data, maxdf = maxdf,
                   mc.cores = cores$m_mc, nthreads = cores$bam_nthreads)
    if(verbose) a = deltime(a, "M-step")

    if(any(lapply(tps, class) == "try-error")) {
      try_errors = try_errors + 1
      if(verbose) {
        cat("\n")
        for(g in 1:k) if(class(tps[[g]])[1] == "try-error")
          print(paste0("Group ", g, " :", attr(tps[[g]], "condition")$message))
        cat("Random reshuffle for next iteration.\n")
      }
      new_group = sample(k, n_id, replace = TRUE)
      changes = sum(new_group != group)
      counts = tabulate(new_group)
    } else {
      ##
      ## E-step:
      ##   predict each id trajectory from each tps model
      pred = mclapply(1:k, pred_g, tps = tps, data = data, mc.cores = cores$e_mc)
      ##   compute loss of each id to each to spline
      loss = mclapply(1:k, mse_g, pred = pred, data = data, mc.cores = cores$e_mc)
      ## classify each id to mse-nearest tps model
      new_group = apply(do.call(cbind, loss), 1, which.min)
      if(verbose) a = deltime(a, " E-step")
      changes = sum(new_group != group)
      counts = tabulate(new_group)
      deviance = sum(unlist(lapply(1:k, function(g) deviance(tps[[g]]))))
      if(verbose)
         cat(" Changes:", changes, "Counts:", counts, "Deviance:", deviance)
    }
    group = new_group
    data$group = as.factor(group[data$id]) # expand group to data frame
    if(changes == 0) break

    ## add outlier treatment as probability of belonging?
    ## use predict.gam(), probably for each group separately? weights??
  }

  if(plot) {
    plot_tps(data, tps)
    if(verbose) a = deltime(a, "\nDone plots")
  }
  if(verbose) deltime(a_0, " trajectories time =")
  list(deviance = deviance, group = group, data_group = data$group, tps = tps,
       iterations = i, try_errors = try_errors, changes = changes)
}

#' Cluster trajectories 
#' 
#' Most users will run the \code{clustra()} function, which takes care of
#' starting values and completes kmeans iteration according to parameters in
#' \code{PL} supplied by \code{clustra_par()} (see \code{demo} directory for
#' examples).
#' 
#' @param data Data frame with response measurements, one per observation.
#' Column names are id, time, response, group.
#' @param k Number of clusters (groups).
#' @param PL A list of lists control structure created by \code{clustra_par}.
#' @param group A vector of initial cluster assignments for unique id's in data.
#' Normally, this is NULL and starts are provided by \code{start_groups()}. 
#' @param verbose Logical to turn on more output during fit iterations.
#' @param plot Logical to indicate 
#' @export
clustra = function(data, k, PL, group = NULL, verbose = FALSE, plot = FALSE) {
  ## get initial group assignments
  if(is.null(group))
    group = start_groups(data, k, nstart = PL$traj_par$starts$ns,
                         nid = PL$traj_par$starts$nid, cores = PL$cores,
                         maxdf = PL$traj_par$maxdf)

  ## provide sequential id's and add initial group assignments to data
  data$id = as.numeric(factor(data$id))
  data$group = group[data$id] # expand to all observations

  ## kmeans iteration to assign id's to groups
  trajectories(data, k, group, PL$traj_par$iter, PL$traj_par$maxdf, plot,
               PL$cores, verbose)
}

