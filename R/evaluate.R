#' allpair_RandIndex is used internally by \code{\link{rand_clustra}} and
#' provides its return value.
#' 
#' Runs \code{\link[MixSim]{RandIndex}} for all pairs of cluster results in its
#' list input and produces a matrix for use by \code{\link{rand_plot}}.
#' 
#' @param results
#' A list with each element packed internally by the
#' \code{\link{rand_clustra}} function.
#' 
#' @return A data frame with \code{\link[MixSim]{RandIndex}} for all pairs from
#' trajectories results. The data frame names and format is intended to be the
#' input for \code{\link{rand_plot}}. Note that all pairs means lower triangle
#' plus diagonal of an all-pairs symmetric matrix.
#' 
#' @export
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
    par(bg = bg, fg = fg, col = fg, col.axis = fg,
        col.lab = fg, col.main = fg, col.sub = fg,
        fig = c(0, ht/wd, 0, 1),
        mar = par()$mar + c(0, 1, 0, -0.5))
    image(x, y, z, zlim, lim, rev(lim), col = var.col, axes = FALSE,
          xlab = "Number of Clusters", ylab = "Number of Clusters",
          main = "Adjusted Rand Index")
    abline(h = (1:K.len) * R.max, lty = 1, col = "black")
    abline(v = (1:K.len) * R.max, lty = 1, col = "black")
    axis(1, at = ((0:(K.len - 1)) + 0.5) * R.max, labels = K.vec)
    axis(2, at = ((0:(K.len - 1)) + 0.5) * R.max, labels = K.vec)
    box()

    par(bg = bg, fg = fg, col = fg, col.axis = fg,
        col.lab = fg, col.main = fg, col.sub = fg, mar = c(5, 0, 4, 2),
        fig = c(ht/wd, 1, 0, .9), new = TRUE)
    z.label = seq(zlim[1], zlim[2], length = n.col)
    image(1, z.label, matrix(z.label, nrow = 1), zlim, c(1, 1), zlim,
          col = var.col, axes = FALSE, xlab = "", ylab = "")
    axis(4)
    box()
  if(is.character(name)) dev.off()
  layout(1)
  invisible(name)
}

#' clustra_sil:
#' Performs \code{\link{clustra}} runs for several k and makes silhouette plots.
#' 
#' Available selection indices are "silhouette" and "rand".
#' Silhouette computes a proxy silhouette index based on distances to cluster
#' centers rather than trajectory pairs. The cost is essentially that of
#' running clustra for several k as this information is
#' available directly from clustra.
#'
#' Rand runs several random start replicates per k, resulting in a much slower
#' algorithm. Then prepares a Rand index comparison between
#' all pairs of clusterings.
#' 
#' @param data
#' The data (see \code{\link{clustra}} description).
#' @param k
#' Vector of k values to try (see \code{.clustra_env}).
#' @param save
#' Logical. When TRUE, save all results as file \code{results.Rdata}.
#' @param verbose
#' Logical. When TRUE, information about each run of clustra is printed.
#' 
#' @return See \code{\link{allpair_RandIndex}}
#' 
#' @export
clustra_sil = function(data, k = clustra_env("ran$ng_vec"),
                         save = FALSE, verbose = FALSE) {
  sil = function(x) {
    ord = order(x)
    ck = ord[1]
    nk = ord[2]
    s = (x[nk] - x[ck])/max(x[ck], x[nk])
    c(ck, nk, s)
  }
  results = vector("list", length(k))
  
  a_rand = deltime()
  for(j in 1:length(k)) {
    kj = k[j]
    a_0 = deltime()
      
    ## Random initial groups to assess stability (NULL gets good starts)
    group = sample(kj, length(unique(data$id)), replace = TRUE)
    f = clustra(data, kj, group, verbose = verbose)
      
    smat = t(apply(f$loss, 1, sil))
    colnames(smat) = c("cluster", "neighbor", "sil_width")
    class(smat) = "silhouette"
    results[[j]] = smat
      
    a = deltime()
    if(verbose) cat(kj, "dev =", f$deviance, "err =", f$try_errors, "ch =",
                    f$changes, "time =", a - a_0, "\n")
  }
  
  ## save object results and parameters
  if(save) save(results, k, file = "clustra_sil.Rdata")
  
  results
}

#' clustra_rand:
#' Performs \code{\link{clustra}} runs for several k and presents a k selection
#' index.
#' 
#' Available selection indices are "silhouette" and "rand".
#' Silhouette computes a proxy silhouette index based on distances to cluster
#' centers rather than trajectory pairs. The cost is essentially that of
#' running clustra for several k as this information is
#' available directly from clustra.
#'
#' Rand runs several random start replicates per k, resulting in a much slower
#' algorithm. Then prepares a Rand index comparison between
#' all pairs of clusterings.
#' 
#' @param data
#' The data (see \code{\link{clustra}} description).
#' @param k
#' Vector of k values to try (see \code{.clustra_env}).
#' @param replicates
#' Number of replicates for each k (see \code{.clustra_env}).
#' @param save
#' Logical. When TRUE, save all results as file \code{results.Rdata}.
#' @param verbose
#' Logical. When TRUE, information about each run of clustra is printed.
#' 
#' @return See \code{\link{allpair_RandIndex}}
#' 
#' @export
clustra_rand = function(data, k = clustra_env("ran$ng_vec"),
                         replicates = clustra_env("ran$replicates"),
                         save = FALSE, verbose = FALSE) {
  results = vector("list", replicates*length(k))
  
  ## Internally, force sequential ids by converting to a factor
  data$id = as.numeric(factor(data$id))
  
  a_rand = deltime()
  for(j in 1:length(k)) {
    kj = k[j]
    for(i in 1:replicates) {
      a_0 = deltime()
      
      ## Random initial groups to assess stability (NULL gets good starts)
      group = sample(kj, length(unique(data$id)), replace = TRUE)
      data$group = group[data$id] # expand group to all data
      f = trajectories(data, kj, group, verbose = verbose)
      if(!is.null( (er = exit_report(f)) )) print(er)
      
      results[[(j - 1)*replicates + i]] = list(k = as.integer(kj), 
                                               rep = as.integer(i),
                                               deviance = f$deviance,
                                               group = f$group)
      
      a = deltime()
      if(verbose) cat(kj, i, "it =", f$iterations, "dev =", f$deviance,
                      "err =", f$try_errors, "ch =", f$changes, "time =",
                      a - a_0, "\n")
    }
  }
  
  ## save object results and parameters
  if(save) save(results, k, replicates, file = "clustra_rand.Rdata")
  
  ## compute and return Rand Index evaluation
  ret = allpair_RandIndex(results)

  ret
}
