#' Function to test information criteria. Not exported and used by internal
#' function `kchoose`.
#' 
#' @param cl 
#' Output from `clustra` function.
#' @param data 
#' A valid data set for `clustra`.
#' @param fn
#' Character file name for output.
#' 
#' @return
#' Numerical value of computed AIC. Also writes a line of computed information 
#' criteria to `fn` file for each k.
#' 
ic_fun = function(cl, data, fn) {
  group = NULL # for data.table R CMD check
  
  k = length(cl$tps)
  N = sum(cl$counts_df)
  
  aic = unlist(lapply(cl$tps, function(x) x$aic))
  deviance = sum(unlist(lapply(cl$tps, function(x) x$deviance)))
  ssq = sum(unlist(lapply(1:k, m = cl$tps, function(i, m) 
     sum((predict(m[[i]], data[group == i]) - data[group == i]$response)^2))))
  df.null = unlist(lapply(cl$tps, function(x) x$df.null))
  df.residual = unlist(lapply(cl$tps, function(x) x$df.residual))
  edf = unlist(lapply(cl$tps, function(x) sum(x$edf)))
  AIC = ssq/N + 2*sum(edf) 
  
  ## (7.30) in ESL (Hastie & Tibshirani)
  sig2 = ssq/(N - sum(edf))
  ## agrees with weighted result from gam/bam of mgcv
  ##sig2w = sum(cl$counts_df * unlist(lapply(cl$tps, function(x) sum(x$sig2))))/N
  AICht = ssq/N + 2*(sum(edf)/N)*sig2
  BIC = ssq/N + log(sum(cl$counts_df))*sum(edf)
  BICht = (N/sig2)*(ssq/N + log(N)*(sum(edf)/N)*sig2) ## (7.36) in ESL (Hastie & Tibshirani)
  BICrams = ssq/sig2 + log(N)*(sum(edf))
  
  ## Oddly, the standard AIC and BIC seem to work well and the adjusted "ht" 
  ## and "rams" do not!? Is the standard correct??
  fmt = "%d %7.0f %6.0f %6.0f %6.0f %6.0f %6.0f %10.0f %10.0f"
  outline = sprintf(fmt, k, sum(aic), AIC, AICht, BIC, BICht, BICrams, deviance,
                    ssq)
  cat(outline, "\n", file = fn, append = TRUE)
  
  AIC
}

#' A test function to evaluate information criteria for several k values. Not
#' exported and only for debugging internal use.
#' 
#' @param K
#' Integer vector of k values to try.
#' @param var
#' A numerical value of noise variance in generated data.
#' @param maxdf 
#' Fitting parameters. See \code{\link{trajectories}}.
#' @param mc 
#' Number of cores to use. Increase up to largest k, or number of cores 
#' available, whichever is less. (On hyperthreaded cores, up to 2x number of 
#' cores.)
#' @param fn 
#' Character file name for output.
#' 
kchoose = function(K, var = 5, maxdf = 10, mc = 1, fn = "ic.txt") {
  library(clustra)
  mc = 4 # If running on a unix or a Mac platform, increase up to 2x # cores
  init = "distant" # maximize median absolute error for initial traj selection
  set.seed(12345)
  data = gen_traj_data(n_id = c( 400, 400, 800, 800, 1600), 
                       types = c(1, 2, 3, 2, 1),
                       intercepts = c(100, 80, 120, 60, 70), m_obs = 25, 
                       s_range = c(-365, -14), e_range = c(0.5*365, 2*365),
                       noise = c(0, var))
  cat("k   aic     AIC  AICht   BIC   BICht   BICrams  deviance ssq\n", file = fn)
  plot_sample(data, layout = c(3,3))
    
  for(k in K) {
    set.seed(1234737)
    cl = clustra(data, k, starts = init, maxdf = maxdf, conv = c(10, 0), 
                 mccores = mc, verbose = TRUE)
    plot_smooths(data, cl$tps)
    ic_fun(cl, data, fn = fn)
  }
  res = data.table::fread(fn)
  
  opar = par(mfrow = c(2, 2))
  on.exit(par(opar))
  plot(AIC ~ k, data = res, type = "l")
  plot(BIC ~ k, data = res, type = "l")
  plot(deviance ~ k, data = res, type = "l")
  plot(ssq ~ k, data = res, type = "l")
}

#' allpair_RandIndex: helper for replicated cluster comparison
#' 
#' Runs \code{\link[MixSim]{RandIndex}} for all pairs of cluster results in its
#' list input and produces a matrix for use by \code{\link{rand_plot}}.
#' Understands replicates within `k` values.
#' 
#' @param results
#' A list with each element packed internally by the
#' \code{\link{clustra_rand}} function with elements: 
#' * `k` - number of clusters
#' * `rep` - replicate number
#' * `deviance` - final deviance
#' * `group` - integer cluster assignments
#' Note that item order is assumed to be the same across all `rep` and `k` but 
#' `group` numbering need not be same. The algorithm only examines if pairs of
#' items are in same or different clusters within each `results` list element.
#' 
#' @return A data frame with \code{\link[MixSim]{RandIndex}} for all pairs from
#' trajectories results. The data frame names and its format is intended to be 
#' the input for \code{\link{rand_plot}}. Note that all pairs is the lower
#' triangle plus diagonal of an all-pairs symmetric matrix.
#' 
allpair_RandIndex = function(results) {
  nr = length(results)
  rand_pairs = vector("list", nr*(nr - 1)/2 + nr)
  irm = 0
  for(i in 1:nr) {
    for(j in i:nr) {
      irm = irm + 1
      irow = results[[i]]
      jrow = results[[j]]
      i.K = as.integer(irow$k)
      i.R = as.integer(irow$rep)
      j.K = as.integer(jrow$k)
      j.R = as.integer(jrow$rep)
      ri = MixSim::RandIndex(irow$group, jrow$group)
      Rand = ri$R
      adjRand = ri$AR
      Fow = ri$F
      Mir = ri$M
      rand_pairs[[irm]] = 
        data.frame(i.K = i.K, i.R = i.R, j.K = j.K, j.R = j.R, 
                   Rand = Rand, adjRand = adjRand, Fow = Fow, Mir = Mir)
    }
  }
  ri_df = do.call(rbind, rand_pairs)
  class(ri_df) = c("RandIndex", class(ri_df))
  ri_df
}


#' clustra_sil: Prepare silhouette plot data for several k or for a previous 
#' clustra run
#' 
#' Performs \code{\link{clustra}} runs for several k and prepares silhouette
#' plot data. Computes a proxy silhouette index based on distances to cluster
#' centers rather than trajectory pairs. The cost is essentially that of
#' running clustra for several k as this information is available directly from
#' clustra. Can also reuse a previous clustra run and produce data for a single
#' silhouette plot.
#' 
#' @param data
#' A data.frame (see the `data` parameter of \code{\link{trajectories}}).
#' Alternatively, the output from a completed `clustra` run can be used, in
#' which case `kv` is left as `NULL`. See Details.
#' @param kv
#' Vector of `clustra` `k` values to run. If `data` is the output from a 
#' completed `clustra` run, leave `kv` as NULL.
#' @param starts
#' See \code{\link{clustra}}.
#' @param mccores
#' See \code{\link{trajectories}}.
#' @param maxdf
#' Fitting parameters. See \code{\link{trajectories}}.
#' @param conv
#' Fitting parameters. See \code{\link{trajectories}}.
#' @param save
#' Logical. When TRUE, save all results as file `clustra_sil.Rdata`.
#' @param verbose
#' Logical. When TRUE, information about each run of clustra is printed.
#' 
#' @details 
#' When given the raw data as the first parameter (input `data` parameter of
#' \code{\link{trajectories}}), `kv` specifies a vector of `k` parameters for 
#' `clustra` and produces data for silhouette plots of each of them. 
#' Alternatively, the input can be the output from a single `clustra` run, in
#' which case data for a single silhouette plot will be made without running 
#' `clustra`.
#' 
#' @return Invisibly returns a list of length `length(kv)`, where each element is
#' a matrix with `nrow(data)` rows and three columns `cluster`, `neighbor`, 
#' `silhouette`. The matrix in each element of this list can be used to draw a 
#' silhouette plot. When the input was a completed `clustra` run, the output is a
#' list with a single element for a single silhouette plot.
#' 
#' @export
clustra_sil = function(data, kv = NULL, starts = "random", mccores = 1, 
                       maxdf = 30, conv = c(10, 0), save = FALSE, 
                       verbose = FALSE) {
  sil = function(x) {
    ord = order(x)
    ck = ord[1]
    nk = ord[2]
    s = (x[nk] - x[ck])/max(x[ck], x[nk])
    c(ck, nk, s)
  }

  ## verify data and kv agreement
  if(is.data.frame(data) && is.null(kv)) {
    cat("clustra_sil: error: must specify kv \n")
    return(NULL)
  } else if(is.matrix(data$loss)) {
      if(is.null(kv)) {
        kv = ncol(data$loss)
      } else if(kv != ncol(data$loss)) {
        cat("clustra: misspecified kv")
        return()
      }
  } else if(is.data.frame(data)) {
    if(length(kv) < 1) {
      cat("clustra: misspecified kv")
      return(NULL)
    }
  } else {
    cat("clustra_sil: error: Expecting data.frame or output from clustra.\n")
    return(NULL)
  }
  
  results = vector("list", length(kv))
  
  for(j in 1:length(kv)) {
    kj = kv[j]
    a_0 = deltime()

    if(is.data.frame(data)) {
      f = clustra(data, kj, starts = starts, mccores = mccores, maxdf = maxdf,
                  conv = conv, verbose = verbose)
    } else {
      f = list(loss = data$loss)
    }

    ## prepare data for silhouette plot
    smat = as.data.frame(t(apply(f$loss, 1, sil)))
    names(smat) = c("cluster", "neighbor", "silhouette")
    smat = smat[order(smat$cluster, -smat$silhouette), ]
    smat$id = factor(rownames(smat), levels = rownames(smat))
    smat$cluster = as.factor(smat$cluster)
    results[[j]] = smat
      
    rm(f, smat)
    gc()
  }
  
  ## save object results and parameters
  if(save) save(results, kv, file = "clustra_sil.Rdata")
  
  invisible(results)
}

#' Function to run trajectories inside mclapply with one core.
#' 
#' @param group
#' Vector of starting group values for unique id's.
#' @param data
#' The data (see \code{\link{clustra}} description).
#' @param k
#' Integer number of clusters.
#' @param maxdf
#' Fitting parameters. See \code{\link{trajectories}}.
#' @param conv
#' Fitting parameters. See \code{\link{trajectories}}.
#' 
#' @return
#' See return of {\code{\link{trajectories}}}.
#' 
traj_rep = function(group, data, k, maxdf, conv) {
  id = ..group = NULL # for data.table R CMD check
  data[, group:=..group[id]] # expand group to all data
  o = trajectories(data, k, group, maxdf, conv = conv, mccores = 1,
                     verbose = FALSE)

  ## reduce output to what's needed downstream (removes big tps and data_group)
  list(deviance = o$deviance, group = o$group, k = o$k, k_cl = o$k_cl,
       iterations = o$iterations, counts = o$counts, 
       counts_df = o$counts_df, changes = o$changes)
}
  
#' clustra_rand: Rand Index cluster evaluation
#' 
#' Performs \code{\link{trajectories}} runs for several *k* and several random
#' start replicates within *k*. Returns a data frame with a Rand Index
#' comparison between all pairs of clusterings. This data frame is typically
#' input to \code{\link{rand_plot}} to produce a heat map with the Adjusted
#' Rand Index results.
#' 
#' @param data
#' The data (see \code{\link{clustra}} description).
#' @param k
#' Vector of k values to try.
#' @param starts
#' See \code{\link{clustra}}.
#' @param mccores
#' Number of cores for replicate parallelism via mclapply.
#' @param replicates
#' Number of replicates for each k.
#' @param maxdf
#' Fitting parameters. See \code{link{trajectories}}.
#' @param conv
#' Fitting parameters. See \code{link{trajectories}}.
#' @param save
#' Logical. When TRUE, save all results as file \code{results.Rdata}.
#' @param verbose
#' Logical. When TRUE, information about each run of clustra (but not iterations
#' within) is printed.
#' 
#' @return See \code{\link{allpair_RandIndex}}.
#' 
#' @export
clustra_rand = function(data, k, starts, mccores, replicates = 10, maxdf = 30,
                        conv = c(10, 0), save = FALSE, verbose = FALSE) {
  id = .GRP = ..group = NULL # for data.table R CMD check
  results = vector("list", replicates*length(k))
  
  ## check for required variables in data
  vnames = c("id", "time", "response")
  if(!is.data.frame(data)) stop("Expecting class data frame and data.table.")
  if(!data.table::is.data.table(data)) data = data.table::as.data.table(data)
  if(!all(vnames %in% names(data))) 
    stop(paste0("Expecting (", paste0(vnames, collapse = ","), ") in data."))
  
  data[, id:=.GRP, by=id] # replace group ids to be sequential
  
  a_rand = deltime()
  n_id = length(unique(data$id))
  for(j in 1:length(k)) {
    kj = k[j]
    
    ## Set replicate starts outside parallel section for reproducibility.
    ## Then run replicates in parallel on mccores. Each replicate is a serial
    ## run of trajectories() with reduced output.
    grp = lapply(rep(kj, replicates), start_groups, 
                 data = data, starts = starts, maxdf = maxdf, conv = conv,
                 mccores = 1, verbose = FALSE)
    tr_lst = parallel::mclapply(grp, traj_rep, data = data, k = kj, 
                                maxdf = maxdf, conv = conv, mc.cores = mccores)
    
    ## generate exit report for each replicate
    xit = lapply(tr_lst, function(x, maxdf, conv)
      if(!is.null( (er = xit_report(x, maxdf, conv)) )) {
        return(er)} else return(NULL), maxdf = maxdf, conv = conv)

    ## package replicates into results for Rand Index computation   
    for(i in 1:replicates) {
      results[[(j - 1)*replicates + i]] = list(k = as.integer(kj), 
                                               rep = as.integer(i),
                                               deviance = tr_lst[[i]]$deviance,
                                               group = tr_lst[[i]]$group)
      if(verbose) 
        cat(kj, i, "iters =", tr_lst[[i]]$iterations, "deviance =", 
            tr_lst[[i]]$deviance, "xit =", xit[[i]], "counts =", 
            tr_lst[[i]]$counts, "changes =", tr_lst[[i]]$changes, "\n")
    }
    rm(tr_lst)
    gc()
  }
  
  ## save object results and parameters
  if(save) save(results, k, replicates, file = "clustra_rand.Rdata")
  
  ## compute and return Rand Index evaluations
  ret = allpair_RandIndex(results)

  ret
}
