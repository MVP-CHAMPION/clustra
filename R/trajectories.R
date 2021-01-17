##
## This is the R translation of trajectories v2.sas code. Possibly, some
## differences exist, which we hope to reconcile in the near future.
##

#' Function to fit thin plate spline (tps) to a group with
#' \code{\link[mgcv]{bam}} from the \code{mgcv} package. Crossvalidation does
#' not know about zero weights, resulting in different smoothing parameters, so
#' subset parameter (rather than zero weights) is used to ensure correct
#' crossvalidation sampling.
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
tps_g = function(data, maxdf, nthreads) {
  tps = mgcv::bam(response ~ s(time, k = maxdf), data = data,
                  discrete = TRUE, nthreads = nthreads)
  tps
}

#' Function to predict for new data based on fitted tps of a group
#'
#' @param tps
#' Output structure of \code{\link[mgcv]{bam}}.
#' @param data
#' See \code{\link{clustra}} description.
#'
pred_g = function(tps, data)
  predict(object = tps, newdata = data, type = "response", se.fit = TRUE)

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
mse_g = function(g, pred, data) { # mean squared error (use unique() to get order)
  ## ctapply assumes id are contiguous within data! can we make sure? 
  data[, esq:=(response - ..pred[[g]]$fit)^2]
#  tt = as.numeric(tapply(data$e_sq, data$id, mean))
  tt = data[, mean(esq), by=id][, 2]
#  ct = as.numeric(iotools::ctapply(esq, data$id, mean))
#  cat(" tt", length(tt), tt[1:10], "\n")
#  cat(" ct", length(ct), ct[1:10], "\n")
  tt
}
#' @rdname mse_g
mxe_g = function(g, pred, data) # maximum error
  as.numeric(tapply(abs(data$response - pred[[g]]$fit), data$id, max))

#' Function to assign starting groups.
#'
#' @param data 
#' Data frame with response measurements, one per observation.
#' Column names are id, time, response, group. Note that id are already
#' sequential starting from 1. This affects expanding group numbers to ids.
#' @param k 
#' Number of clusters (groups).
#' @param fp
#' Fitting parameters. See \code{link{trajectories}}.
#' @param cores 
#' 
#' @param verbose 
#' Turn on more output for debugging.
#'
#' @importFrom methods is
#' @export
start_groups = function(data, k, fp,
                        cores = c(e_mc = 1, m_mc = 1, nthreads = 1, blas = 1),
                        verbose = FALSE) {
  if(verbose) cat("\nStarts : ")

  ## Prepare test data for diversity criterion
  test_data = data.frame(id = rep(0, k*2*fp$maxdf),
                time = rep(seq(data[, min(time)], data[, max(time)],
                               length.out = 2*fp$maxdf), times = k),
                response = rep(NA, k*2*fp$maxdf),
                group = rep(1:k, each = 2*fp$maxdf))

  ## Samples k*idperstart id's. Picks tps fit with best deviance after one
  ##   iteration among random starts. Choosing from samples increases diversity
  ##   of fits (sum of distances between group fits).
  ## Then classifies all ids based on fit from best sample
  max_div = 0
  start_id = sort(sample(data[, unique(id)], k*fp$idperstart)) # for speed and diversity!
  data_start = data[id %in% start_id]

  data_start[, id:=.GRP, by=id] # replace id's to be sequential
  ## Note: since id are not sequential, am I losing correspondence to ids? No,
  ##  because spline models do not know about ids. The classification process
  ##  only needs who has same id, not what is the id.)

  ## evaluate starts on id sample
  for(i in 1:fp$starts) {
    group = sample(k, k*fp$idperstart, replace = TRUE) # random groups
    data_start[, group:=..group[id]] # expand group to all data

    fp$iter = 1  # single iter for starts only here
    f = trajectories(data_start, k, group, fp, cores = cores, verbose = verbose)
    er = xit_report(f, fp)
    if(verbose && !is.null((er))) cat(" ", er)
    
    diversity = sum(dist(
        do.call(rbind, lapply(parallel::mclapply(f$tps, pred_g,
                                       data = test_data, mc.cores = cores["e_mc"]),
                              function(x) as.numeric(x$fit)))
    ))
    if(diversity > max_div) {
      max_div = diversity
      best_tps_cov = f$tps
    }
    if(verbose) cat(" Diversity:", round(diversity, 2), "")
  }

  if(max_div == 0) return(NULL) # all starts failed!

  ## expand best start sample fit to full set of id's
  pred = parallel::mclapply(best_tps_cov, pred_g, data = data, 
                            mc.cores = cores["e_mc"])
  ## compute mse for each id
  loss = parallel::mclapply(1:k, mse_g, pred = pred, data = data,
                  mc.cores = cores["e_mc"]) # !!! loss is not data dimensioned!!!!
  ## classify id to group with min mse
  group = apply(do.call(cbind, loss), 1, which.min)
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
#' Group numbers corresponding to sequential ids.
#' @param fp
#' List with addtional fitting parameters: *maxdf* basis dimension of smooth
#' term (see \code{\link[mgcv]{s}}, parameter k, in package \code{mgcv}),
#' *iter* maximum number of EM iterations, *starts* number of trial random starts
#' from which the b  est is chosen, *idperstart* number of `id`s to sample for
#' evaluating random starts, *retry_max* number of restarts if iteration
#' encounters a cluster too small for \code{\link[mgcv]{bam}} fitting with the
#' given *maxdf*.
#' @param cores
#' Vector of length 4. Core allocations to various sections of code: `e_mc` dores to use
#' by `mclapply` of e_mc (expectation over k), `m_mc` cores to use by 
#' `mclapply` of m_mc (maximization across k), `nthreads` threads to use by
#' \link[mgcv]{bam}, `blas` cores to use by OpenBLAS.
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
  
  ## make sure that data is a data.table
  if(!data.table::is.data.table(data)) data = data.table::as.data.table(data)
  
  if(verbose) a = a_0 = deltime(a)
  if(data[, min(id)] != 1 | data[, max(id)] != length(group)){
    cat("\ntrajectories: Expecting sequential id's starting from 1.\n")
  }

  ## get number of unique id
  n_id = data.table::uniqueN(data[, id])

  ## EM algorithm to cluster ids into k groups
  ## iterates fitting a thin plate spline (tps) center to each group (M-step)
  ##      regroups each id to nearest tps center (E-step)
  for(i in 1:fp$iter) {
    if(verbose) cat("\n", i, "")
    ##
    ## M-step:
    ##   Estimate tps model parameters for each cluster from id's in that cluster
    datg = parallel::mclapply(1:k, function(g) data[group == g])
    tps = parallel::mclapply(datg, tps_g, maxdf = fp$maxdf,
                             mc.cores = cores["m_mc"], nthreads = cores["nthreads"])
    if(verbose) a = deltime(a, "M-step")

    ## break if any bam fit fails
    if(!all(sapply(tps, is, class = "bam"))) {
      changes = counts_df = counts = loss = NULL
      break
    }

    ##
    ## E-step:
    ##   predict each id trajectory by each model
    pred = parallel::mclapply(tps, pred_g, data = data, mc.cores = cores["e_mc"])
    ##   compute loss of each id to each model
    loss = parallel::mclapply(1:k, mse_g, pred = pred, data = data,
                              mc.cores = cores["e_mc"])
    loss = do.call(cbind, loss)

    ## classify each id to model with smallest loss
    new_group = apply(loss, 1, which.min)
    if(verbose) a = deltime(a, " E-step")

    ## evaluate results and update groups
    changes = sum(new_group != group)
    counts = tabulate(new_group, k)
    deviance = sum(unlist(lapply(1:k, function(g) deviance(tps[[g]]))))
    if(verbose)
       cat(" Changes:", changes, "Counts:", counts, "Deviance:", deviance)
    group = new_group
    data[, group:=..new_group[id]] # expand group to data frame
    counts_df = data[, .N, by=group][, N]
    
    ## break if converged or if any cluster gets too small for fp$maxdf
    if(changes == 0 || any(counts_df <= fp$maxdf) || length(counts_df) < k) break
  }

  if(verbose) deltime(a_0, "\n trajectories time =")
  list(deviance = deviance, group = group, loss = loss, data_group = data[, group],
       tps = tps, iterations = i, counts = counts, counts_df = counts_df,
       changes = changes)
}

#' xit_report examines trajectories output to name what was concluded, such as
#' convergence, maximum iterations reached, a zero cluster, etc.
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
  if(!is.null(cl$counts) && any(cl$counts == 0))
    xit = c(xit, "zerocluster")
  if(!is.null(cl$changes) && cl$changes == 0)
    xit = c(xit, "converged")
  if(!all(sapply(cl$tps, is, class = "bam")))
    xit = c(xit, "bamfail")
  if(cl$iterations >= fp$iter && fp$iter != 1)
    xit = c(xit, "max iter")
  xit
}

#' Cluster trajectories
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
#' @param group
#' A vector of initial cluster assignments for unique id's in data.
#' Normally, this is NULL and good starts are provided by
#' \code{link{start_groups}}.
#' @param fp
#' Fitting parameters. See \code{link{trajectories}}.
#' @param verbose
#' Logical to turn on more output during fit iterations.
#' 
#' @examples
#' set.seed(13)
#' data = gen_traj_data(n_id = 200, m_obs = 20, s_range = c(-50, -10),
#'               e_range = c(3*365, 10*365), reference = 100)
#' cl = clustra(data, k = 3, verbose = TRUE)
#' names(cl)
#' cl$tps
#'
#' @export
clustra = function(data, k, group = NULL,
                   fp = list(maxdf = 30, iter = 8, starts = 3, idperstart = 20,
                   retry_max = 3), cores = c(e_mc = 1, m_mc = 1, nthreads = 1,
                                             blas = 1), verbose = FALSE) {
  ## check for required variables in data
  vnames = c("id", "time", "response")
  if(!is.data.frame(data)) stop("Expecting a data frame.")
  ## further, make sure that data is a data.table (convert if not)
  if(!data.table::is.data.table(data)) data = data.table::as.data.table(data)
  if(!all(vnames %in% names(data))) 
    stop(paste0("Expecting (", paste0(vnames, collapse = ","), ") in data."))
  
  data[, id:=.GRP, by=id] # replace id's to be sequential
  if(!is.null(group)) data[, group:=..group[id]] # expand group to all data

  for(retry in 1:fp$retry_max) {
    ## get initial group assignments
    while(is.null(group)) {
      group = start_groups(data, k, fp, verbose = verbose)
      data[, group:=..group[id]] # expand group to all data
      if(any(data[, .(low = .N < ..fp$maxdf, by=group)][, low])) {
        group = NULL
        cat("\nRepeating starts due to cluster observations under fp$maxdf ...")
      }
    }

    ## kmeans iteration to assign id's to groups
    cl = trajectories(data, k, group, fp, cores, verbose = verbose)
    er = xit_report(cl, fp)
    if(verbose && !is.null(er)) cat(" ", er)
    if( any( c("converged", "max iter") %in% er ) )break
    cat("\n Restarting clustra. Error exit. \n")
  }
  cl$retry = retry

  cl
}
