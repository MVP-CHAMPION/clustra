##
## This is the R translation of trajectories v2.sas code
##

## Function to fit thin plate spline (tps) to a group with gam from the mgcv
## package. Crossvalidation does not know about zero weights, resulting in
## different smoothing parameters, so subset parameter is used to ensure 
## correct crossvalidation sampling.
tps_g = function(g, data, maxdf, nthreads) {
  tps = mgcv::bam(response ~ s(time, k = maxdf), data = data,
                  subset = data$group == g, discrete = TRUE, nthreads = nthreads)
  if(class(tps)[1] == "try-error") print(attr(tps, "condition")$message)
  tps
}

#' Function to predict for new data based on fitted tps of a group
#' 
#' @param tps
#' Output structure of [mgcv::bam()]
#' @param data 
#' See [clustra()] description.
#' 
pred_g = function(tps, data)
  predict(object = tps, newdata = data, type = "response", se.fit = TRUE)
## Loss functions for each id to a group tps center
mse_g = function(g, pred, data) # mean squared error
  as.numeric(tapply((data$response - pred[[g]]$fit)^2, data$id, mean))
mxe_g = function(g, pred, data) # maximum error
  as.numeric(tapply(abs(data$response - pred[[g]]$fit), data$id, max))

#' Function to assign starting groups.
#' 
#' @param data Data frame with response measurements, one per observation.
#' Column names are id, time, response, group. Note that id are already
#' sequential starting from 1. This affects expanding group numbers to ids.
#' @param k Number of clusters (groups).
#' "bestrand" is implemented.
#' @param nstart (see clu component of .clustra_env)
#' @param nid (see clu component of .clustra_env)
#' @param maxdf (see clu component of .clustra_env)
#' @param cores (see cor component of .clustra_env)
#' @param verbose Turn on more output for debugging.
#' 
#' @importFrom methods is
#' @export
start_groups = function(data, k,
                        nstart = clustra_env("clu$starts"),
                        nid = clustra_env("clu$idperstart"),
                        maxdf = clustra_env("clu$maxdf"),
                        cores = clustra_env("cor"),
                        verbose = FALSE) {
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
  data_start = data[match(data$id, start_id, nomatch = 0) > 0, ]

  data_start$id = as.numeric(factor(data_start$id)) 
  ## Note: since id are not sequential (Am I losing correspondence to ids? No,
  ##  because spline modesl do not know about ids. The classification process 
  ##  only needs who has same id, not what is the id.)

  ## evaluate starts on id sample
  for(i in 1:nstart) {
    group = sample(k, k*nid, replace = TRUE) # random groups
    data_start$group = group[data_start$id] # expand group to all responses

    f = trajectories(data_start, k, group, iter = 1,  # single iter for starts!
                     verbose = verbose)
    if(!all(sapply(f$tps, is, class = "bam"))) next

    diversity = sum(dist(
        do.call(rbind, lapply(parallel::mclapply(f$tps, pred_g,
                                       data = test_data, mc.cores = cores$e_mc),
                              function(x) as.numeric(x$fit)))
    ))
    if(diversity > max_div) {
      max_div = diversity
      best_tps_cov = f$tps
    }
    if(verbose) cat(round(diversity, 2), "")
  }

  if(max_div == 0) return(NULL) # all starts failed!
    
  ## expand best start sample fit to full set of id's
  pred = parallel::mclapply(best_tps_cov, pred_g, data = data, mc.cores = cores$e_mc)
  ## compute mse for each id
  loss = parallel::mclapply(1:k, mse_g, pred = pred, data = data,
                  mc.cores = cores$e_mc) # !!! loss is not data dimensioned!!!!
  ## classify id to group with min mse
  group = apply(do.call(cbind, loss), 1, which.min)
  if(verbose) cat("\nStart counts", tabulate(group, k), "->", max_div, "")
  group ## report group for each id
}


#' Cluster longitudinal trajectories of a response variable.
#' 
#' Trajectory means are splines fit to all ids in a cluster. Typically, this
#' function is called by \code{clustra()}.
#' 
#' @param data 
#' Data frame with response measurements, one per observation. Column
#' names are \code{id}, \code{time}, \code{response}, \code{group}. Note that
#' id must be already sequential starting from 1. This affects expanding group 
#' numbers to ids.
#' @param k 
#' Number of clusters (groups)
#' @param group 
#' Group numbers corresponding to sequential ids.
#' @param iter 
#' Maximum iterations in mgcv::bam (see .clustra_env)
#' @param maxdf 
#' Maximum degrees of freedom for tps in mgcv::bam (see .clustra_env)
#' @param cores 
#' List with cores allocation to various sections (see .clustra_env)
#' @param verbose 
#' Logical, whether to produce debug output.
#' 
#' @importFrom stats predict
#' @importFrom methods is
#' 
#' @author George Ostrouchov and David Gagnon
#' 
#' @export
trajectories = function(data, k, group,
                        iter = clustra_env("clu$iter"),
                        maxdf = clustra_env("clu$maxdf"),
                        cores = clustra_env("cor"),
                        verbose = FALSE) {
  
  if(verbose) a = a_0 = deltime(a)
  openblasctl::openblas_set_num_threads(cores$blas)
  if(max(data$id) != length(group)){
    cat("\ntrajectories: group id's may be corrupted or empty!\n")
    browser()
  }

  ## get number of unique id
  n_id = length(unique(data$id))

  ## EM algorithm to cluster ids into k groups
  ## iterate fitting a thin plate spline center to each group (M-step)
  ##         regroup each id to nearest tps center (E-step)
    exit = "max iter"
  for(i in 1:iter) {
    if(verbose) cat("\n", i, "")
    ##
    ## M-step:
    ##   Estimate tps model parameters for each cluster from id's in that cluster
    tps = parallel::mclapply(1:k, tps_g, data = data, maxdf = maxdf,
                   mc.cores = cores$m_mc, nthreads = cores$bam_nthreads)
    if(verbose) a = deltime(a, "M-step")

    if(!all(sapply(tps, is, class = "bam"))) {
      exit = "try-error"
      changes = -1
      break
    } else {
      ##
      ## E-step:
      ##   predict each id trajectory from each tps model
      pred = parallel::mclapply(tps, pred_g, data = data, mc.cores = cores$e_mc)
      ##   compute loss of each id to each to spline
      loss = parallel::mclapply(1:k, mse_g, pred = pred, data = data,
                                mc.cores = cores$e_mc)
      
      ## classify each id to mse-nearest tps model
      new_group = apply(do.call(cbind, loss), 1, which.min)
      if(verbose) a = deltime(a, " E-step")
      changes = sum(new_group != group)
      counts = tabulate(new_group, k)
      deviance = sum(unlist(lapply(1:k, function(g) deviance(tps[[g]]))))
      if(verbose)
         cat(" Changes:", changes, "Counts:", counts, "Deviance:", deviance)
      group = new_group
      data$group = as.factor(group[data$id]) # expand group to data frame
    }
    if(changes == 0) {
      exit = "converged"
      break
    }
    if(any(tabulate(data$group, k) < maxdf)) {
      exit = "failed: not enough points in cluster"
      break
    }
  }

  if(verbose) deltime(a_0, "\n trajectories time =")
  if(verbose) cat("\n exit:", exit, "\n")
  list(deviance = deviance, group = group, data_group = data$group, tps = tps,
       iterations = i, changes = changes, exit = exit)
}

#' Cluster trajectories 
#' 
#' Most users will run the \code{clustra()} function, which takes care of
#' starting values and completes kmeans iteration according to parameters in
#' \code{.clustra_env} environment (See vignette("clustra_basic.Rmd") for
#' more details).
#' 
#' @param data 
#' Data frame with response measurements, one per observation.
#' Column names are id, time, response, group.
#' @param k 
#' Number of clusters (groups).
#' @param group 
#' A vector of initial cluster assignments for unique id's in data.
#' Normally, this is NULL and starts are provided by \code{start_groups()}. 
#' @param verbose 
#' Logical to turn on more output during fit iterations.
#' @param rngkind
#' Character string giving random number generator type (see [RNGkind()]).
#' @param seed
#' Seed for generating random starting initial cluster assignments.
#' 
#' @details In addition to the shown parameters, detailed clustering and core
#' allocation parameters are in .clustra_env environment that can be controlled
#' with \code{clustra_env} function.
#' 
#' @export
clustra = function(data, k, group = NULL, verbose = FALSE,
                   rngkind = clustra_env("clu$rngkind"),
                   seed = clustra_env("clu$seed")) {
  ## set RNG and its seed
  rng_prev = RNGkind(rngkind)
  set.seed(seed)

  for(retry in 1:clustra_env("clu$retry_max")) {
    ## get initial group assignments
    while(is.null(group)) {
      group = start_groups(data, k, verbose = verbose)
      ## provide sequential id's and add initial group assignments to data
      data$id = as.numeric(factor(data$id))
      data$group = group[data$id] # expand to all observations
      if(any(tabulate(data$group, k) < clustra_env("clu$maxdf"))) {
        group = NULL
        cat("\nRepeating starts due to cluster observations under clu$maxdf ...")
      }
    }
    
    ## kmeans iteration to assign id's to groups
    cl = trajectories(data, k, group, verbose = verbose)
    if(cl$exit == "converged" | cl$exit == "max iter") break
    cat("\n Restarting clustra due to", cl$exit, ".\n")
  }
  cl$retry = retry
  
  ## restore RNG
  RNGkind(rng_prev[1])
  
  cl
}
