##
## This is the R translation of trajectories v2.sas code. Possibly, some
## differences exist, which we hope to reconcile in the near future.
##

#' Function to fit thin plate spline (tps) to a group with
#' \code{\link[mgcv]{bam}} from the \code{mgcv} package. 
#' Crossvalidation does not know about zero weights, resulting in
#' different smoothing parameters, so subset parameter (rather than
#' zero weights) is used to ensure correct crossvalidation sampling.
#'
#' @param g
#' Integer group number.
#' @param data
#' See \code{\link{clustra}} description.
#' @param maxdf
#' See \code{\link{trajectories}} description.
#' @param nthreads
#' Controls \code{\link[mgcv]{bam}} threads.
#'
tps_g = function(g, data, maxdf, nthreads) {
  if(nrow(data[[g]]) > 0) {
    return(mgcv::bam(response ~ s(time, k = maxdf), data = data[[g]],
                  discrete = TRUE, nthreads = nthreads))
  } else {
    return(NULL)
  }
}

#' Function to predict for new data based on fitted tps of a group
#'
#' @param tps
#' Output structure of \code{\link[mgcv]{bam}}.
#' @param data
#' See \code{\link{clustra}} description.
#'
pred_g = function(tps, data) {
  if(is.null(tps)) {
    return(NULL)
  } else {
    predict(object = tps, newdata = data, type = "response", se.fit = TRUE)
  }
}

#' Loss functions
#'
#' \code{mse_g} Computes mean-squared loss for each group.
#' \code{mxe_g} is maximum loss.
#'
#' @param g
#' Integer group number.
#' @param pred
#' List. Element g is output of \code{\link{pred_g}} for group g.
#' @param data
#' See \code{\link{trajectories}}
#'
mse_g = function(g, pred, data) { # mean squared error
  esq = response = ..pred = id = NULL # for data.table R CMD check
  if(is.null(pred[[g]])) {
    return(NULL)
  } else {
    data[, esq:=(response - ..pred[[g]]$fit)^2]
    tt = as.numeric(unlist(data[, mean(esq), by=id][, 2]))
    return(tt)
  }
}
#' @rdname mse_g
mxe_g = function(g, pred, data) # maximum error
  as.numeric(tapply(abs(data$response - pred[[g]]$fit), data$id, max))

#' Function to check if non-zero group members have enough data for spline fit
#' degrees of freedom.
#' 
#' @param group
#' Group membership vector for each id.
#' @param loss
#' A matrix with rows of computed loss values of each id across all models
#' @param data
#' A data.table with data. See \code{\link{trajectories}}.
#' @param fp
#' Fitting parameters. See \code{\link{trajectories}}.
#' 
#' @return 
#' Returns the vector of group membership of id's either unchanged or changed
#' to have sufficient data in non-zero groups.
#' 
check_df = function(group, loss, data, fp) {
  ..group = id = NULL # for data.table R CMD check
  counts_df = data[, tabulate(group)]
  while(any(counts_df > 0 & counts_df <= fp$maxdf)) {  # need to reallocate
    low_gp = which(counts_df > 0 & counts_df <= fp$maxdf)[1]
    loss[, low_gp] = rep(Inf, nrow(loss)) # set to Inf to zero group
    move_id = (group == low_gp)
    group[move_id] = apply(loss[move_id, , drop = FALSE], 1, which.min)
    data[, group:=..group[id]]  # push group change into data
    counts_df = data[, tabulate(group)]
  }
  group
}

#' Function to assign starting groups.
#' 
#' If only one start, a random assignment is done. If more than one start, 
#' picks tps fit with best deviance after one iteration among random starts. 
#' Choosing from samples increases diversity of fits (sum of distances between 
#' group fits). Then classifies all ids based on fit from best sample.
#'
#' @param data 
#' Data frame with response measurements, one per observation.
#' Column names are id, time, response, group. Note that id are already
#' sequential starting from 1. This affects expanding group numbers to ids.
#' @param k 
#' Number of clusters (groups).
#' @param starts
#' A vector of length two, giving the number of start values to try and the
#' number of ids per cluster to evaluate the starts (If the number of ids is 
#' less than 1, use all data and do not subset data for starts.). If more than 
#' one start, the most diverse after one trajectories iteration on a sample of 
#' data is used. The default is *c(1, 0)*, meaning that one random start is 
#' used with all the data.
#' @param fp
#' Fitting parameters. See \code{\link{trajectories}}.
#' @param cores
#' A vector of core assignments for multicore components. See
#' @param verbose 
#' Turn on more output for debugging.
#'
#' @importFrom methods is
#' @export
start_groups = function(data, k, starts, fp,
                        cores = c(e_mc = 1, m_mc = 1, nthreads = 1, blas = 1),
                        verbose = FALSE) {
  time = id = .GRP = ..group = NULL # for data.table R CMD check
  
  if(verbose) cat("\nStarts : ")

  ## Prepare test data for diversity criterion
  test_data = data.frame(id = rep(0, k*2*fp$maxdf),
                time = rep(seq(data[, min(time)], data[, max(time)],
                               length.out = 2*fp$maxdf), times = k),
                response = rep(NA, k*2*fp$maxdf),
                group = rep(1:k, each = 2*fp$maxdf))

  ## Sample a subset of the data (for speed and increased diversity of starts)
  max_div = 0
  if(starts[2] > 0) {
    start_id = sort(sample(data[, unique(id)], k*starts[2]))
    data_start = data[id %in% start_id]
    data_start[, id:=.GRP, by=id] # replace id's to be sequential in the subset
  } else {
    data_start = data
  }
  n_id = data_start[, data.table::uniqueN(id)]
  
  ## evaluate starts on the subset
  for(i in 1:starts[1]) {
    group = sample(k, n_id, replace = TRUE) # random groups
    data_start[, group:=..group[id]] # replicate group numbers to ids

    fp$iter = 1  # use single fit iteration for starts
    f = trajectories(data_start, k, group, fp, cores = cores, verbose = verbose)
    er = xit_report(f, fp)
    if(verbose && !is.null((er))) cat(" ", er)
    
    diversity = sum(dist( # compute diversity of fit across test data
        do.call(rbind, lapply(parallel::mclapply(f$tps, pred_g,
                                       data = test_data, mc.cores = cores["e_mc"]),
                              function(x) as.numeric(x$fit)))
    ))
    if(diversity > max_div) { # record the biggest diversity across starts
      max_div = diversity
      best_tps_cov = f$tps
    }
    if(verbose) cat(" Diversity:", round(diversity, 2), "")
  }

  if(max_div == 0) return(NULL) # all starts failed!

  ## predict from best sample fit to full set of id's
  pred = parallel::mclapply(best_tps_cov, pred_g, data = data, 
                            mc.cores = cores["e_mc"])
  ## compute loss for each id wrt model of each cluster
  loss = parallel::mclapply(1:k, mse_g, pred = pred, data = data,
                  mc.cores = cores["e_mc"])
  loss = do.call(cbind, loss) # combine list into matrix columns
  group = apply(loss, 1, which.min) # set group as closest cluster mean
  data[, group:=..group[id]] # replicate group numbers within ids
  group = check_df(group, loss, data, fp) # check sufficient df
  
  if(verbose) cat("\nStart counts", tabulate(group, k), "->", max_div, "")
  group ## report group for each id
}


#' Cluster longitudinal trajectories of a response variable.
#'
#' Trajectory means are splines fit to all ids in a cluster. Typically, this
#' function is called by \code{\link{clustra}}.
#'
#' @param data
#' Data frame with response measurements, one per observation. Column
#' names are \code{id}, \code{time}, \code{response}, \code{group}. Note that
#' id must be already sequential starting from 1. This affects expanding group
#' numbers to ids.
#' @param k
#' Number of clusters (groups)
#' @param group
#' Vector of initial group numbers corresponding to ids.
#' @param fp
#' List with additional fitting parameters: *maxdf* basis dimension of smooth
#' term (see \code{\link[mgcv]{s}}, parameter k, in package \code{mgcv}),
#' *iter* maximum number of EM iterations.
#' @param cores
#' Vector of length 4. Core allocations to various sections of code: *e_mc*
#' cores to use by *mclapply* of e_mc (expectation over k), *m_mc* cores to use
#' by *mclapply* of m_mc (maximization across k), *nthreads* threads to use by
#' \link[mgcv]{bam}, cores to use by OpenBLAS.
#' @param verbose
#' Logical, whether to produce debug output.
#'
#' @importFrom stats predict
#' @importFrom methods is
#'
#' @author George Ostrouchov and David Gagnon
#'
#' @export
trajectories = function(data, k, group, fp,
                        cores = c(e_mc = 1, m_mc = 1, nthreads = 1, blas = 1),
                        verbose = FALSE) {
  if(verbose) a = a_0 = deltime(a)

  id = ..new_group = ..group = NULL # for data.table R CMD check
  
  ## make sure that data is a data.table
  if(!data.table::is.data.table(data)) data = data.table::as.data.table(data)
  k_cl = k  # start with k clusters

  ## get number of unique id
  n_id = data[, data.table::uniqueN(id)]

  ## EM algorithm to cluster ids into k groups
  ## iterates fitting a thin plate spline (tps) center to each group (M-step)
  ##      regroups each id to nearest tps center (E-step)
  for(i in 1:fp$iter) {
    if(verbose) cat("\n", i, "")
    ##
    ## M-step: Estimate model parameters for each cluster
    ##   
    counts = tabulate(group)
    datg = parallel::mclapply(1:k_cl, function(g) data[group == g])
    nz = which(sapply(datg, nrow) > 0) # nonzero groups
    k_cl = length(nz) # reset number of clusters to nonzeros only
    tps = parallel::mclapply(nz, tps_g, data = datg, maxdf = fp$maxdf,
                             mc.cores = cores["m_mc"], 
                             nthreads = cores["nthreads"])
    if(verbose) a = deltime(a, "M-step")

    ##
    ## E-step:
    ##   predict each id's trajectory with each model
    pred = parallel::mclapply(tps, pred_g, data = data, 
                              mc.cores = cores["e_mc"])
    ##   compute loss of all id's to all groups (models)
    ##   TODO better parallel balance by skipping empty groups
    loss = parallel::mclapply(1:k_cl, mse_g, pred = pred, data = data,
                              mc.cores = cores["e_mc"])
    loss = do.call(cbind, loss) # NULL loss elements go away (removes 0 groups)
    ## classify each id to model with smallest loss (Expected group)
    ## Note: Groups are renumbered as NULL loss is compressed in above do.call.
    new_group = apply(loss, 1, which.min) # already without original zero groups
    if(verbose) a = deltime(a, " E-step")

    ## evaluate results and update groups
    changes = sum(new_group != group) # may exaggerate due to renumbering
    counts = tabulate(new_group)
    deviance = sum(unlist(lapply(1:k_cl, function(g) deviance(tps[[g]]))))
    if(verbose)
       cat(" Changes:", changes, "Counts:", counts, "Deviance:", deviance)
    group = new_group
    data[, group:=..new_group[id]] # replicate group to ids
    group = check_df(group, loss, data, fp)
    data[, group:=..group[id]]
    counts_df = data[, tabulate(group)]

    ## break if converged
    if(changes == 0) break
  }

  if(verbose) deltime(a_0, "\n trajectories time =")
  list(deviance = deviance, group = group, loss = loss, k = k, k_cl = k_cl,
       data_group = data[, group], tps = tps, iterations = i, counts = counts,
       counts_df = counts_df, changes = changes)
}

#' xit_report examines trajectories output to name what was concluded, such as
#' convergence, maximum iterations reached, a zero cluster, etc. Multiple
#' conclusions are possible as not all are mutually exclusive.
#' 
#' @param cl
#' Output structure from trajectories() function
#' @param fp
#' Fitting parameters. See \code{link{trajectories}}.
#' 
#' @export
xit_report = function(cl, fp) {
  xit = NULL
  if(!is.null(cl$counts_df) && any(cl$counts_df < fp$maxdf))
    xit = c(xit, "undermaxdf")
  if(cl$k_cl < cl$k)
    xit = c(xit, "zerocluster")
  if(!is.null(cl$changes) && cl$changes == 0)
    xit = c(xit, "converged")
  if(cl$iterations >= fp$iter && fp$iter != 1)
    xit = c(xit, "max-iter")
  xit
}

#' Cluster trajectoriesß
#'
#' Most users will run the wrapper \code{\link{clustra}} function, which takes
#' care of starting values. See vignette("clustra_basic.Rmd") for
#' more details.
#'
#' @param data
#' Data frame or, preferably, also a data.table with response measurements, one
#' response per observation. Required variables are (id, time, response).
#' Other variables are ignored.
#' @param k
#' Number of clusters
#' @param starts
#' A vector of length two. See \code{\link{start_groups}}.
#' @param group
#' A vector of initial cluster assignments for unique id's in data.
#' Normally, this is NULL and good starts are provided by
#' \code{\link{start_groups}}.
#' @param fp
#' Fitting parameters. See \code{\link{trajectories}}.
#' @param cores
#' A vector of core assignments for multicore components. See
#' \code{\link{trajectories}}. 
#' @param verbose
#' Logical to turn on more output during fit iterations.
#' 
#' @examples
#' set.seed(13)
#' data = gen_traj_data(n_id = 100, m_obs = 20, s_range = c(-50, -10),
#'               e_range = c(2*365, 5*365))
#' fp = list(maxdf = 20, iter = 5, starts = 2, idperstart = 20, retry_max = 3)
#' cl = clustra(data, k = 3, fp = fp, verbose = TRUE)
#' names(cl)
#' cl$tps
#'
#' @export
clustra = function(data, k, starts = c(1, 0), group = NULL,
                   fp = list(maxdf = 30, iter = 10),
                   cores = c(e_mc = 1, m_mc = 1, nthreads = 1, blas = 1),
                   verbose = FALSE) {
  id = .GRP = ..group = NULL # for data.table R CMD check
  
  ## check for required variables in data
  vnames = c("id", "time", "response")
  if(!is.data.frame(data)) stop("Expecting a data frame.")
  if(!data.table::is.data.table(data)) data = data.table::as.data.table(data)
  if(!all(vnames %in% names(data))) 
    stop(paste0("Expecting (", paste0(vnames, collapse = ","), ") in data."))
  
  data[, id:=.GRP, by=id] # replace id's to be sequential
  if(!is.null(group)) data[, group:=..group[id]] # replicate group to ids

  ## get initial group assignments
  group = start_groups(data, k, starts, fp, verbose = verbose)

  ## kmeans iteration to assign id's to groups
  cl = trajectories(data, k, group, fp, cores, verbose = verbose)
  er = xit_report(cl, fp)
  if(verbose && !is.null(er)) cat(" ", er)

  cl
}
