#' Various functions to evaluate results of clustering, like RandIndex
#' computation and plots of spline fits.
#'
#' Plots spline fit (modify to use acutal bam results), possibly with data
#' TODO use actual tps fit rather than redo it with stat_smooth
#' 
#' @section Details:
#'
#' @param dat Data frame with variables *time*, *response*, *color*, *group*
#' @param points If true, plot all the data points.
plot_tps = function(dat, points = TRUE) {
  p = ggplot(dat, aes(time, response, color = group, group = group))
  if(points) p = p + geom_point(alpha = 0.1)
  p = p + stat_smooth(method = "gam", formula = y ~ s(x, k = maxdf), size = 1)
  p = p + theme_bw()
  print(p)
}

#' Use RandIndex to evaluate number of clusters
#' @param results List with results from trajectories() function
#' @return A data frame with RandIndex for all pairs from trajectories results.
#' Note all pairs means lower triangle plus diagonal of an all-pairs symmetric
#' matrix.
allpair_RandIndex = function(results) {
  library(MixSim)
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
      ri = RandIndex(irow$group, jrow$group)
      Rand = ri$R
      adjRand = ri$AR
      Fow = ri$F
      Mir = ri$M
      rand_pairs[[irm]] = 
        data.frame(i.K = i.K, i.R = i.R, j.K = j.K, j.R = j.R, 
                   Rand = Rand, adjRand = adjRand, Fow = Fow, Mir = Mir)
    }
  }
  do.call(rbind, rand_pairs)
}

#' Rand index matrix plot from Technometrics paper
#' Sorts replicates within cluster K
#' Assumes K starts from 2
#' Author: Wei-Chen Chen (wcsnow@gmail.com)
## 
#' @param rand_pairs A data frame with columns of cluster assignments
rand_plot = function(rand_pairs, name) {
  library(RColorBrewer)

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
  var.col = colorRampPalette(brewer.pal(9, "YlOrRd"))(n.col)

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
  pdf(name, width = wd, height = ht, bg = bg)
    par(bg = bg, fg = fg, col = fg, col.axis = fg,
        col.lab = fg, col.main = fg, col.sub = fg,
        fig = c(0, ht/wd, 0, 1),
        mar = par()$mar + c(0, 1, 0, -0.5))
    image(x, y, z, zlim, lim, rev(lim), col = var.col, axes = FALSE,
          xlab = "Number of Clusters", ylab = "Number of Clusters",
          main = "Adjusted Rand Index")
    abline(h = (1:K.len) * R.max, lty = 1, col = "black")
    abline(v = (1:K.len) * R.max, lty = 1, col = "black")
    axis(1, at = ((0:(K.len - 1)) + 0.5) * R.max, label = K.vec)
    axis(2, at = ((0:(K.len - 1)) + 0.5) * R.max, label = K.vec)
    box()

    par(bg = bg, fg = fg, col = fg, col.axis = fg,
        col.lab = fg, col.main = fg, col.sub = fg, mar = c(5, 0, 4, 2),
        fig = c(ht/wd, 1, 0, .9), new = TRUE)
    z.label = seq(zlim[1], zlim[2], length = n.col)
    image(1, z.label, matrix(z.label, nrow = 1), zlim, c(1, 1), zlim,
          col = var.col, axes = FALSE, xlab = "", ylab = "")
    axis(4)
    box()
  dev.off()
  layout(1)
}

#' Performs clustra runs for several k and several random start replicates per k.
#' Then prepares a Rand index comparison between all pairs of clusterings, that
#' are displayed in a matrix plot.
#' @param data The data (see clustra description).
#' @param k Vector of k values to run.
#' @param PL A list data structure (a list of lists) giving all needed
#' parameters for the clustra runs.
#' @param save Logical. When TRUE, save all results as file results.Rdata
#' @param verbose Logical. When TRUE, information about each run of clustra is
#' printed.
#' @export
rand_clustra = function(data, k, PL, save = FALSE, verbose = FALSE) {
  set.seed(PL$traj_par$seed)
  replicates = PL$traj_par$replicates
  
  a_rand = deltime()
  results = vector("list", replicates*length(k_vec))
  for(j in 1:length(k_vec)) {
    kj = k[j]
    for(i in 1:replicates) {
      a_0 = deltime()
      
      f = clustra(data, kj, PL)
      results[[(j - 1)*replicates + i]] = list(k = as.integer(kj), 
                                               rep = as.integer(i),
                                               deviance = f$deviance,
                                               group = f$group)
      data$group = as.factor(f$group[data$id])
      a = deltime()
      if(verbose) cat(kj, i, "it =", f$iterations, "dev =", f$deviance,
                      "err =", f$try_errors, "ch =", f$changes, "time =",
                      a - a_0, "\n")
    }
  }
  
  ## save object results and parameters
  if(save) save(results, k, maxdf, starts, iter, replicates,
                file = "results.Rdata")
  
  ## plot Rand Index evaluation
  RandIndex_pairs = allpair_RandIndex(results)
  rand_plot(RandIndex_pairs, name = "adjRand_mat.pdf")
  if(verbose) a = deltime(a_rand, "\nTotal rand_cluster time")
}

