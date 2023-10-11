# Make sure data.table knows we know we're using it
.datatable.aware = TRUE

#' Fits a thin plate spline to a single group with
#' \code{\link[mgcv]{bam}}. 
#' 
#' @description 
#' Fits a thin plate spline to a single group (one list element) in data with
#' \code{\link[mgcv]{bam}}. Uses data from only one group rather than a
#' zero weights approach. Zero weights would result in
#' incorrect crossvalidation sampling.
#'
#' @param g
#' Integer group number.
#' @param data
#' A list of group-separated data using lapply with 
#' \code{data.table::copy(data[group == g])} from original data in
#' \code{\link{clustra}} description.
#' @param maxdf
#' See \code{\link{trajectories}} description.
#' @param nthreads
#' Controls \code{\link[mgcv]{bam}} threads.
#' @return 
#' Returns an object of class "gam". See \code{\link[mgcv]{bam}} value. 
#' If group data has zero rows, NULL is returned instead.
tps_g = function(g, data, maxdf, nthreads) {
  if(nrow(data[[g]]) > 0) {
    return(mgcv::bam(response ~ s(time, k = maxdf), data = data[[g]],
                  discrete = TRUE, nthreads = nthreads))
  } else {
    return(NULL)
  }
}

#' Function to predict for new data based on fitted gam object.
#'
#' @param tps
#' Output structure of \code{\link[mgcv]{bam}}.
#' @param newdata
#' See \code{\link{clustra}} description of data parameter.
#' @return 
#' A numeric vector of predicted values corresponding to rows of newdata. 
#' If gam object is NULL, NULL is returned instead.
pred_g = function(tps, newdata) {
  if(is.null(tps)) {
    return(NULL)
  } else {
    return(as.vector(mgcv::predict.bam(object = tps, newdata = newdata, type = "response",
                             newdata.guaranteed = TRUE)))
  }
}

#' Various Loss functions used internally by clustra
#'
#' @param pred
#' Vector of predicted values.
#' @param id
#' Integer vector of group assignments.
#' @param response
#' Vector of response values.
#' @return 
#' A numeric value. For mse_g(), returns the mean-squared error. For mae_g(),
#' returns mean absolute error. For mme_g(), returns median absolute error.
#' For mxe_g(), returns the
#' maximum absolute error.
mse_g = function(pred, id, response) {
  if(is.null(pred)) {
    return(NULL)
  } else {
    esq = (response - pred)^2
    DT = data.table::data.table(esq, id)
    tt = as.numeric(unlist(DT[, mean(esq), by=id][, 2]))
    return(tt)
  }
}
#' @rdname mse_g
mae_g = function(pred, id, response) {
  if(is.null(pred)) {
    return(NULL)
  } else {
    esq = abs(response - pred)
    DT = data.table::data.table(esq, id)
    tt = as.numeric(unlist(DT[, mean(esq), by=id][, 2]))
    return(tt)
  }
}
#' @rdname mse_g
mme_g = function(pred, id, response) {
  if(is.null(pred)) {
    return(NULL)
  } else {
    esq = abs(response - pred)
    DT = data.table::data.table(esq, id)
    tt = as.numeric(unlist(DT[, median(esq), by=id][, 2]))
    return(tt)
  }
}
#' @rdname mse_g
mxe_g = function(pred, id, response) { # maximum error
  if(is.null(pred)) {
    return(NULL)
  } else {
    eabs = abs(response - pred)
    DT = data.table::data.table(eabs, id)
    tt = as.numeric(unlist(DT[, max(eabs), by=id][, 2]))
    return(tt)
  }
}


#' Checks if non-empty groups have enough data for spline fit
#' degrees of freedom.
#' 
#' @param group
#' An integer vector of group membership for each id.
#' @param loss
#' A matrix with rows of computed loss values of each id across all models
#' as columns.
#' @param data
#' A data.table with data. See \code{\link{trajectories}}.
#' @param maxdf
#' Fitting parameters. See \code{\link{trajectories}}.
#' @details 
#' When a group has insufficient data for \code{maxdf}, its nearest model loss
#' values are set to \code{Inf}, and new nearest model is assigned. The check
#' repeats until all groups have sufficient data.
#' @return 
#' Returns the vector of group membership of id's either unchanged or changed
#' to have sufficient data in non-zero groups.
#' 
check_df = function(group, loss, data, maxdf) {
  ..group = id = NULL # for data.table R CMD check
  counts_df = data[, tabulate(group)]
  while(any(counts_df > 0 & counts_df <= maxdf)) {  # need to reallocate
    low_gp = which(counts_df > 0 & counts_df <= maxdf)[1]
    loss[, low_gp] = rep(Inf, nrow(loss)) # set to Inf to zero-out group
    move_id = (group == low_gp)
    group[move_id] = apply(loss[move_id, , drop = FALSE], 1, which.min)
    data[, group:=..group[id]]  # push group change into data
    counts_df = data[, tabulate(group)]
  }
  return(group)
}

#' Function to assign starting groups.
#' 
#' Either a random assignment of k approximately equal size clusters or a
#' FastMap-like algorithm that sequentially selects k distant ids
#' from those that have more than the median number of observations. TPS fits 
#' to these ids are used as cluster centers for a starting group assignment.
#' A user supplied starting assignment is also possible.
#'
#' @param data 
#' Data.table with response measurements, one per observation.
#' Column names are id, time, response, group. Note that \code{id}s are assumed
#' sequential starting from 1. This affects expanding group numbers to ids.
#' @param k 
#' Number of clusters (groups).
#' @param starts
#' Type of start groups generated. See \code{\link{clustra}}.
#' @param maxdf
#' Fitting parameters. See \code{\link{trajectories}}.
#' @param conv
#' Fitting parameters. See \code{\link{trajectories}}.
#' @param mccores
#' See \code{\link{trajectories}}.
#' @param verbose 
#' Turn on more output for debugging. Values 0, 1, 2, 3 add more output. 2 and
#' 3 produce graphs during iterations - use carefully!
#' @return 
#' An integer vector corresponding to unique `id`s, giving group number
#' assignments.
#' 
#' For `distant`, each sequential selection takes an id that has the largest
#' minimum distance from smooth TPS fits (<= 5 deg) of previous selections. 
#' The distance of an id to a single TPS is the median absolute error across 
#' the id time points. Distance of an id to a set of TPS is the minimum of 
#' the individual distances. We pick the id that has the maximum of such 
#' a minimum of medians.
#'
#' @importFrom methods is
#' @importFrom stats median
start_groups = function(k, data, starts, maxdf, conv, mccores = 1,
                        verbose = FALSE) {
  time = response = id = .GRP = ..group = NULL # for data.table R CMD check
  if(verbose) a = a_0 = deltime(a)
  
  ## get number of unique id
  n_id = data[, data.table::uniqueN(id)]
  
  if(starts == "distant") {
    ## Use a FastMap-like algorithm to generate a starting configuration
    ## among id's with above-median number of observations
    sid = vector("integer", k)
    sfit = vector("list", k)
    sdata = vector("list", k)
    tid = tabulate(data$id)
    if(length(tid) != n_id) print("start_groups: Wrong # ids?!")
    mtid = median(tid)
    gmid = which(tid > mtid) # only ids with more than median observations
    gdata = data[id %in% gmid]
  
    ## Start random. No need to remove first random.
    sid[1] = sample(gmid, 1)
    sdata[[1]] = gdata[id==sid[1]]
    sfit[[1]] = mgcv::gam(response ~ s(time), data = sdata[[1]], maxdf = 5)
    pred = pred_g(sfit[[1]], gdata)
    loss = mme_g(pred, gdata$id, gdata$response) # median absolute error

    ## Replace with farthest from first and remove farthest from available
    imax = which.max(loss) # largest median absolute error
    sid[1] = gmid[imax]
    gmid = gmid[-imax]
    sdata[[1]] = gdata[id==sid[1]]
    sfit[[1]] = mgcv::gam(response ~ s(time), data = sdata[[1]], maxdf = 5)
    gdata = gdata[id!=sid[1]]

    for(i in 2:k) {
      loss = rep(Inf, length(gmid)) 
      for(j in 1:(i - 1)) { # a vector of minimum distance to previous i - 1
        pred = pred_g(sfit[[j]], gdata)
        loss = pmin(loss, mme_g(pred, gdata$id, gdata$response)) 
      }
      imax = which.max(loss) # id with largest min dist to previous (maximin)
      sid[i] = gmid[imax] # set next
      gmid = gmid[-imax] # remove from available ids
      sdata[[i]] = gdata[id==sid[i]] # get its values and fit spline next
      sfit[[i]] = mgcv::gam(response ~ s(time), data = sdata[[i]], maxdf = 5)
      gdata = gdata[id!=sid[i]] # remove its data
    }
    
    ## plot selected ids (undocumented debug verbose)
    if(verbose > 2) plot_smooths(data = data, fits = sfit, select.data = sdata)
     
    ## now classify all data to the distant subjects
    newdata = force(as.data.frame(data))
    pred = parallel::mclapply(sfit, pred_g, newdata = newdata, mc.cores = mccores)
    rm(newdata, gdata, sdata, sfit)
    gc()

    ## compute loss for each id wrt model of each cluster
    loss = parallel::mclapply(pred, mse_g, id = force(data[, id]),
                              response = force(data[, response]),
                              mc.cores = mccores)
    loss = do.call(cbind, loss) # combine list into matrix columns
    group = apply(loss, 1, which.min) # set group as closest cluster mean
    
    if(verbose) cat("\n distant ids: ", sid, "")
  } else if(starts == "random") {
    group = sample(k, n_id, replace = TRUE)
  } else if(is.integer(starts) && 
            length(starts) == n_id && 
            all(sort(unique(starts)) == 1:k)) {
    group = starts
  } else {
    cat("starts invalid. Proceeding with random groups.\n")
    group = sample(k, n_id, replace = TRUE)
  }

  data = data[, group:=..group[id]] # replicate group numbers within ids
  
  # Check if sufficient observations. Drop and reassign group if not.
  group = check_df(group, loss, data, maxdf)
    
  if(verbose) {
    cat("  Start counts", tabulate(group, k), "")
    deltime(a_0, "  Starts time: ")
  }
  return(group) ## return group for each id
}

#' Cluster longitudinal trajectories over time.
#'
#' Performs k-means clustering on continuous `response` measured over `time`, 
#' where each mean is defined by a thin plate spline fit to all points in a
#' cluster. Typically, this function is called by \code{\link{clustra}}.
#'
#' @param data
#' Data table or data frame with response measurements, one per observation.
#' Column names are `id`, `time`, `response`, `group`. Note that
#' `id`s must be sequential starting from 1. This affects expanding group
#' numbers to `id`s.
#' @param k
#' Number of clusters (groups)
#' @param group
#' Vector of initial group numbers corresponding to `id`s.
#' @param maxdf
#' Integer. Basis dimension of smooth term. See \code{\link[mgcv]{s}} function
#' parameter `k`, in package `mgcv`.
#' @param conv
#' A vector of length two, `c(iter, minchange)`, where `iter` is the maximum
#' number of EM iterations and `minchange` is the minimum percentage of subjects
#' changing group to continue iterations. Setting `minchange` to zero continues
#' iterations until no more changes occur or `maxiter` is reached.
#' @param mccores
#' Integer number of cores to use by `mclapply` sections. Parallelization is 
#' over `k`, the number of clusters.
#' @param verbose
#' Logical, whether to produce debug output. A value > 1 will plot tps fit lines
#' in each iteration.
#' @param ...
#' See \code{\link{clustra}} for allowed `...` parameters.
#' 
#' @return 
#' A list with components
#' * `deviance` - The final deviance in each cluster added across clusters.
#' * `group` - Integer vector of group assignments corresponding to unique `id`s.
#' * `loss` - Numeric matrix with rows corresponding to unique `id`s and one 
#' column for each cluster. Each entry is the mean squared loss for the data in
#' the `id` relative to the cluster model.
#' * `k` - An integer giving the requested number of clusters.
#' * `k_cl` - An integer giving the converged number of clusters. Can be 
#' smaller than `k` when some clusters become too small for degrees of freedom
#' during convergence. 
#' * `data_group` - An integer vector, giving group assignment as expanded into
#' all `id` time points.
#' * `tps` - A list with `k_cl` elements, each an object returned by the 
#' `mgcv::bam` fit of a cluster thin plate spline model.
#' * `iterations` - An integer giving the number of iterations taken.
#' * `counts` - An integer vector giving the number of `id`s in each cluster.
#' * `counts_df` - An integer vector giving the total number of observations in
#' each cluster (sum of the number of observations for `id`s belonging to the 
#' cluster).
#' * `changes` - An integer, giving the number of `id`s that changed clusters in
#' the last iteration. This is zero if converged.
#' 
#' @importFrom stats predict
#' @importFrom methods is
#'
#' @author George Ostrouchov and David Gagnon
#'
#' @export
trajectories = function(data, k, group, maxdf, conv = c(10, 0), mccores = 1,
                        verbose = FALSE, ...) {
  if(verbose) a = a_0 = deltime(a)
  xargs = list(...)

  time = response = id = ..new_group = ..group = NULL # for data.table R CMD check
  
  ## make sure that data is a data.table
  if(!data.table::is.data.table(data)) data = data.table::as.data.table(data)
  k_cl = k  # start with k clusters

  ## get number of unique id
  n_id = data[, data.table::uniqueN(id)]

  ## EM algorithm to cluster ids into k groups
  ## iterates fitting a thin plate spline (tps) center to each group (M-step)
  ##      regroups each id to nearest tps center (E-step)
  for(i in 1:conv[1]) {
    if(verbose) cat("\n", i, "")
#    gc() # Tighten up memory use in each
    ##
    ## M-step: Estimate model parameters for each cluster ---------------------
    ##   
    if(verbose) cat("(M-")
    datg = parallel::mclapply(1:k_cl, 
                              function(g) data.table::copy(data[group == g]),
                              mc.cores = 1)
    if(verbose && any(sapply(datg, is.null))) cat("*C*")
    if(verbose) cat("1")
    nz = which(sapply(datg, nrow) > 0) # nonzero groups
    k_cl = length(nz) # reset number of clusters to nonzeros only
    if(verbose) cat("2")
    tps = parallel::mclapply(nz, tps_g, data = datg, maxdf = maxdf,
                             mc.cores = 1, nthreads = 1)
    if(verbose && any(sapply(tps, is.null))) cat("*F*")
    if(verbose) a = deltime(a, "3)")

    ##
    ## E-step: predict (classify) each id to a model --------------------------
    ## 
    if(verbose) cat(" (E-")
    newdata = force(as.data.frame(data[, list(time, response)]))
    pred = parallel::mclapply(tps, pred_g, newdata = newdata, 
                              mc.cores = mccores)
    if(verbose && any(sapply(pred, is.null))) cat("*P*")
    if(verbose) cat("1")
    rm(newdata)
#    gc() # tighten up memory before next mclapply
    ##   compute loss of all id's to all groups (models)
    ##   TODO better parallel balance by skipping empty groups
    if(verbose) cat("2")
    loss = parallel::mclapply(pred, mse_g,
                             id = force(data[, id]),
                             response = force(data[, response]),
                             mc.cores = mccores)
    rm(pred)
    if(verbose && any(sapply(loss, is.null))) cat("*E*")
    if(verbose) cat("3")
    loss = do.call(cbind, loss) # NULL loss elements go away (removes 0 groups)
    if(verbose) cat("4")
    ## classify each id to model with smallest loss (Expected group)
    ## Note: Groups are renumbered as NULL loss is compressed in above do.call.
    new_group = apply(loss, 1, which.min) # already without original zero groups
    if(verbose) a = deltime(a, "5)")

    ##
    ## evaluate results and update groups -------------------------------------
    ## 
    changes = sum(new_group != group) # may exaggerate due to renumbering
    counts = tabulate(new_group)
    deviance = sum(unlist(lapply(1:k_cl, function(g) deviance(tps[[g]]))))
    if(verbose)
       cat(" Changes:", changes, "Counts:", counts, "Deviance:", deviance)
    group = new_group
    data[, group:=..new_group[id]] # replicate group to ids
    group = check_df(group, loss, data, maxdf)
    data[, group:=..group[id]]
    counts_df = data[, tabulate(group)]
    
    if(verbose > 1) plot_smooths(data, tps, max.data = 0, ...)

    ## break if converged
    if(100*changes/sum(counts) <= conv[2]) break
  }

  ## Compute AIC and BIC
  N = nrow(data)
  ssq = unlist(lapply(1:length(tps), m = tps, function(i, m)
    sum((predict(m[[i]], data[group == i]) - data[group == i]$response)^2)))
  edf = round(sum(unlist(lapply(tps, function(x) sum(x$edf)))), 2)
  AIC = round(sum(ssq)/N + 2*edf, 2) 
  BIC = round(sum(ssq)/N + log(N)*edf, 2)

  msg = paste0("\n AIC:", AIC, " BIC:", BIC, ' edf:', edf, "  Total time: ")
  if(verbose) deltime(a_0, msg)
  return(list(deviance = deviance, group = group, loss = loss, k = k, k_cl = k_cl,
       data_group = data[, group], tps = tps, iterations = i, counts = counts,
       counts_df = counts_df, changes = changes, AIC = AIC, BIC = BIC))
}

#' xit_report 
#' 
#' Examines trajectories output to name what was concluded, such as
#' convergence, maximum iterations reached, a zero cluster, etc. Multiple
#' conclusions are possible as not all are mutually exclusive.
#' 
#' @param cl
#' Output structure from \code{\link{trajectories}} function
#' @param maxdf
#' Fitting parameters. See \code{\link{trajectories}}.
#' @param conv
#' Fitting parameters. See \code{\link{trajectories}}.
#' @return 
#' NULL or a character vector of exit criteria satisfied.
#' 
xit_report = function(cl, maxdf, conv) {
  xit = NULL
  if(!is.null(cl$counts_df) && any(cl$counts_df < maxdf))
    xit = c(xit, "undermaxdf")
  if(cl$k_cl < cl$k)
    xit = c(xit, "zerocluster")
  if(!is.null(cl$changes) && 100*cl$changes/sum(cl$counts) <= conv[2])
    xit = c(xit, "converged")
  if(cl$iterations >= conv[1] && conv[1] != 1)
    xit = c(xit, "max-iter")
  return(xit)
}

#' Cluster longitudinal trajectories over time
#'
#' The usual top level function for clustering longitudinal trajectories. After
#' initial setup, it calls \code{\link{trajectories}} to perform k-means 
#' clustering on continuous `response` measured over `time`, where each mean 
#' is defined by a thin plate spline fit to all points in a cluster. See
#' `clustra_vignette.Rmd` for examples of use.
#'
#' @param data
#' Data frame or, preferably, also a data.table with response measurements, one
#' response per observation. Required variables are (id, time, response).
#' Other variables are ignored.
#' @param k
#' Number of clusters
#' @param starts
#' One of c("random", "distant") or an integer vector 
#' with values 1:k corresponding to unique ids of starting cluster assignments.
#' For "random", starting clusters are assigned
#' at random. 
#' For "distant", a FastMap-like algorithm selects k distant ids to
#' which TPS models are fit and used as starting cluster centers to which ids 
#' are classified. Only id with more than median number of time points are
#' used. Distance from an id to a TPS model is median absolute difference
#' at id time points. Starting with a random id, distant ids are selected
#' sequentially as the id with the largest minimum absolute distance to 
#' previous selections (a maximin concept). The first random selection is 
#' discarded and the next k selected ids are kept. Their TPS fits become the 
#' first cluster centers to 
#' which all ids are classified. See comments in code and
#' DOI: 10.1109/TPAMI.2005.164 for the FastMap analogy.
#' @param maxdf
#' Fitting parameters. See \code{\link{trajectories}}.
#' @param conv
#' Fitting parameters. See \code{\link{trajectories}}.
#' @param mccores
#' See \code{\link{trajectories}}. 
#' @param verbose
#' Logical to turn on more output during fit iterations.
#' @param ...
#' Additional parameters of optional plotting under `verbose = 2`. At this time,
#' only `xlim` and `ylim` are allowed.
#' 
#' @return 
#' A list returned by \code{\link{trajectories}} plus one more element `ido`,
#' giving the original id numbers is invisibly returned. Invisible returns are
#' useful for repeated runs that explore verbose clustra output.
#' 
#' @examples
#' set.seed(13)
#' data = gen_traj_data(n_id = c(50, 100), types = c(1, 2), 
#'                      intercepts = c(100, 80), m_obs = 20, 
#'                      s_range = c(-365, -14), e_range = c(0.5*365, 2*365))
#' cl = clustra(data, k = 2, maxdf = 20, conv = c(5, 0), verbose = TRUE)
#' tabulate(data$group)
#' tabulate(data$true_group)
#'
#' @export
clustra = function(data, k, starts = "random", maxdf = 30, conv = c(10, 0),
                   mccores = 1, verbose = FALSE, ...) {
  id = .GRP = .SD = ..group = NULL # for data.table R CMD check

  ## check for required variables in data
  xargs = list(...)
  unkn_param = !(names(xargs) %in% c("ylim", "xlim"))
  if(any(unkn_param))
    stop(paste("unknown parameter(s) ...:", 
               paste(names(xargs)[unkn_param], collapse = " ")))
  vnames = c("id", "time", "response")
  if(!is.data.frame(data)) stop("Expecting a data frame.")
  if(!data.table::is.data.table(data)) data = data.table::as.data.table(data)
  data.table::setalloccol(data, n = 1) # room for group column in start_groups()
  
  if(!all(vnames %in% names(data))) 
    stop(paste0("Expecting (", paste0(vnames, collapse = ","), ") in data."))
  
  ## replace id's to be sequential for faster cluster to data expansions
  ##   keep old id numbers in `ido`
  ido = data[, .SD[1], id][, id]
  data[, id:=.GRP, by=id]
  
  n_id = data[, data.table::uniqueN(id)]
  
  ## Get initial group assignments. Populate into data in-place via data.table
  group = start_groups(k, data, starts, maxdf, conv, mccores,
                       verbose = verbose)
  
  ## Perform k-means iteration for groups
  cl = trajectories(data, k, group, maxdf, conv, mccores, verbose = verbose, ...)
  er = xit_report(cl, maxdf, conv)
  if(verbose && !is.null(er)) cat(" ", er, "\n")
  cl$ido = ido

  invisible(cl)
}
