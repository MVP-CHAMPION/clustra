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

#' Matrix plot of Rand Index comparison of replicated clusters
#' 
#' @param rand_pairs
#' A data frame result of \code{\link{allpair_RandIndex}}
#' @param name 
#' Character string file name for pdf plot. If omitted or NULL, plot will 
#' render to current graphics device.
#'
#' @return
#' Invisible. Full path name of file with plot.
#'
#' @references
#' Wei-chen Chen, George Ostrouchov, David Pugmire, Prabhat, and Michael Wehner.
#' 2013. A Parallel EM Algorithm for Model-Based Clustering Applied to the
#' Exploration of Large Spatio-Temporal Data. Technometrics, 55:4, 513-523.
#'
#' Sorts replicates within cluster K
#' Assumes K starts from 2
#' 
#' @author
#' Wei-Chen Chen and George Ostrouchov
#' 
#' @importFrom 
#' grDevices colorRampPalette dev.off pdf
#' @importFrom 
#' graphics abline axis box image layout par
#' @export
rand_plot = function(rand_pairs, name = NULL) {
  K.vec = unique(unlist(rand_pairs[, c("i.K", "j.K")]))
  K.max = max(K.vec)
  K.len = length(K.vec)
  R.vec = unique(unlist(rand_pairs[, c("i.R", "j.R")]))
  R.len = length(R.vec)
  R.max = max(R.vec)
  n.col = K.len*R.len

  ## axes
  x = 0:((K.max - 1) * R.max)
  y = x
  
  ## fill square array
  z = array(NA, c(max(x), max(x)))
  diag(z) = 1
  var.z = rand_pairs$adjRand
  for(i in 1:nrow(rand_pairs)){ # allows non-sequential K values
    Kbase = function(val) (which(K.vec == val) - 1)*R.max
    z.i = Kbase(rand_pairs$i.K[i]) + rand_pairs$i.R[i]
    z.j = Kbase(rand_pairs$j.K[i]) + rand_pairs$j.R[i]
    z[z.i, z.j] = z[z.j, z.i] = var.z[i]
  }
  lim = range(x)
  zlim = range(c(rand_pairs$adjRand, 1))
  ## made by colors = brewer.pal(9, "YlOrRd") # removed dependency
  colors = c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A",
             "#E31A1C", "#BD0026", "#800026")
  var.col = colorRampPalette(colors)(n.col)
  
  ## Reorder within K for vis effect
  for(i in 1:K.len) { 
    iv = (i - 1) * R.max + (1:R.max)
    for(j in R.len:1) { # using later rows for breaking ties in earlier rows
      ic = (i - 1) * R.max + j
      id.z = order(z[iv, ic], decreasing = TRUE)
      tmp.z = z[iv, ]
      z[iv, ] = tmp.z[id.z,]
      tmp.z = z[, iv]
      z[, iv] = tmp.z[, id.z]
    }
  }

  bg = "transparent"
  fg = "black"
  wd = 6.5
  ht = 6
  if(is.character(name)) pdf(name, width = wd, height = ht, bg = bg)
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar)) # restore old par when exit
    par(bg = bg, fg = fg, col = fg, col.axis = fg, col.lab = fg, col.main = fg,
        col.sub = fg, fig = c(0, ht/wd, 0, 1),
        mar = par()$mar + c(0, 1, 0, -0.5))
    image(x, y, z, zlim, lim, rev(lim), col = var.col, axes = FALSE,
          xlab = "Number of Clusters", ylab = "Number of Clusters",
          main = "Adjusted Rand Index")
    abline(h = (1:K.len) * R.max, lty = 1, col = "black")
    abline(v = (1:K.len) * R.max, lty = 1, col = "black")
    axis(1, at = ((0:(K.len - 1)) + 0.5) * R.max, labels = K.vec)
    axis(2, at = ((0:(K.len - 1)) + 0.5) * R.max, labels = K.vec)
    box()

    par(bg = bg, fg = fg, col = fg, col.axis = fg, col.lab = fg, col.main = fg,
        col.sub = fg, mar = c(5, 0, 4, 2), fig = c(ht/wd, 1, 0, .9), new = TRUE)
    z.label = seq(zlim[1], zlim[2], length = n.col)
    image(1, z.label, matrix(z.label, nrow = 1), zlim, c(1, 1), zlim,
          col = var.col, axes = FALSE, xlab = "", ylab = "")
    axis(4)
    box()
  if(is.character(name)) dev.off()
  layout(1)
  invisible(name)
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
clustra_sil = function(data, kv = NULL, mccores = 1, maxdf = 30, conv = c(10, 0),
                       save = FALSE, verbose = FALSE) {

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
      f = clustra(data, kj, mccores = mccores, maxdf = maxdf, conv = conv,
                verbose = verbose)
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
  id = ..group = NULL
  data[, group:=..group[id]] # expand group to all data
  trajectories(data, k, group, maxdf, conv, 1, verbose = FALSE)
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
#' Logical. When TRUE, information about each run of clustra is printed.
#' 
#' @return See \code{\link{allpair_RandIndex}}.
#' 
#' @export
clustra_rand = function(data, k, mccores, replicates = 10, maxdf = 30,
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
    
    ## Set starting groups outside parallel section to guarantee reproducibility
    grp = lapply(rep(kj, replicates), sample.int, size = n_id, replace = TRUE)
    f = parallel::mclapply(grp, traj_rep, data = data, k = kj, maxdf = maxdf,
                           conv = conv, mc.cores = mccores)
    fer = lapply(f, function(f, maxdf, conv)
      if(!is.null( (er = xit_report(f, maxdf, conv)) )) {
        return(er)} else return(NULL), maxdf = maxdf, conv = conv)
    for(i in 1:replicates) {
      results[[(j - 1)*replicates + i]] = list(k = as.integer(kj), 
                                               rep = as.integer(i),
                                               deviance = f[[i]]$deviance,
                                               group = f[[i]]$group)
      if(verbose) 
        cat(kj, i, "iters =", f[[i]]$iterations, "deviance =", f[[i]]$deviance,
            "xit =", fer[[i]], "counts =", f[[i]]$counts, "changes =",
            f[[i]]$changes, "\n")
    }
    rm(f)
    gc()
  }
  
  ## save object results and parameters
  if(save) save(results, k, replicates, file = "clustra_rand.Rdata")
  
  ## compute and return Rand Index evaluation
  ret = allpair_RandIndex(results)

  ret
}
