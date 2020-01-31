## Various functions to evaluate results of clustering, like RandIndex
## computation and plots of spline fits.
##
#' Plots spline fit (modify to use acutal bam results), possibly with data
#' 
#' @section Details:
#'
#' @param dat  
plot_tps = function(dat, data = TRUE) {
  p = ggplot(dat, aes(time, response, color = group, group = group)) +
  if(data) p = p + geom_point(alpha = 0.1)
  p = p + stat_smooth(method = "gam", formula = y ~ s(x, k = maxdf), size = 1)
  p = p + theme_bw()
  print(p)
}

#' Use RandIndex to evaluate number of clusters
#' @param results List with results from trajectories() function
#' @value A data frame with RandIndex for all pairs from trajectories results.
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
      i.K = as.integer(irow$ng)
      i.R = as.integer(irow$rep)
      j.K = as.integer(jrow$ng)
      j.R = as.integer(jrow$rep)
      ri = RandIndex(irow$group, jrow$group)
      Rand = ri$R
      adjRand = ri$AR
      Fow = ri$F
      Mir = ri$M
      rand_pairs[[irm]] = data.frame(i.K = i.K, i.R = i.R, j.K = j.K, j.R = j.R, 
                                   Rand = Rand, adjRand = adjRand, Fow = Fow, Mir = Mir)
    }
  }
  do.call(rbind, rand_pairs)
}

## Rand index matrix plot from Technometrics paper
## Sorts replicates within cluster K
## Assumes K starts from 2
## Author: Wei-Chen Chen
## 
#' @param rand_pairs A data frame with columns 
rand_plot = function(rand_pairs) {
  library(RColorBrewer)

  n.col = max(rand_pairs$i.K)*max(rand_pairs$i.R)
  K.vec = unique(unlist(rand_pairs[, c("i.K", "j.K")]))
  K.max = max(K.vec)
  K.len = length(K.vec)
  unique.R = unique(unlist(rand_pairs[, c("i.R", "j.R")]))
  R.vec = max(unique.R)
  R.max = max(R.vec)

  x = 1:((K.max - 1) * R.vec)
  y = x
  z = array(NA, c(max(x), max(x)))
  diag(z) = 1
  var.z = rand_pairs$adjRand
  for(i in 1:nrow(rand_pairs)){
    z.i = (rand_pairs$i.K[i] - 2) * R.vec + rand_pairs$i.R[i]
    z.j = (rand_pairs$j.K[i] - 2) * R.vec + rand_pairs$j.R[i]
    z[z.i, z.j] = z[z.j, z.i] = var.z[i]
  }
  lim = range(x)
  zlim = range(c(rand_pairs$adjRand, 1))
  var.col = colorRampPalette(brewer.pal(9, "YlOrRd"))(n.col)

  ## Reorder within K for vis effect
  for(i in 1:K.len) { 
    id.z = order(z[(i - 1) * R.max + (1:R.max), (i - 1) * R.max + 1], decreasing = TRUE)
    tmp.z = z[(i - 1) * R.max + (1:R.max), ]
    z[(i - 1) * R.max + (1:R.max), ] = tmp.z[id.z,]
    tmp.z = z[, (i - 1) * R.max + (1:R.max)]
    z[, (i - 1) * R.max + (1:R.max)] = tmp.z[, id.z]
  }

  bg = "transparent"
  fg = "black"

  pdf("all_sort1.pdf", width = 6, height = 6, bg = bg)
    par(bg = bg, fg = fg, col = fg, col.axis = fg,
        col.lab = fg, col.main = fg, col.sub = fg)
    image(x, y, z, zlim, lim, rev(lim), col = var.col, axes = FALSE,
          xlab = "Number of Clusters", ylab = "Number of Clusters",
          main = "Adjusted Rand Index")
    abline(h = (1:(K.max - 1)) * R.vec + 0.5, lty = 1, col = "black")
    abline(v = (1:(K.max - 1)) * R.vec + 0.5, lty = 1, col = "black")
    axis(1, at = ((0:(K.max - 2)) + 0.5) * R.vec, label = 2:K.max)
    axis(2, at = ((0:(K.max - 2)) + 0.5) * R.vec, label = 2:K.max)
    box()
  dev.off()

  pdf("all_sort2.pdf", width = 0.8, height = 5, bg = bg)
    par(bg = bg, fg = fg, col = fg, col.axis = fg,
        col.lab = fg, col.main = fg, col.sub = fg, mar = c(5, 1, 4, 2))
    z.label = seq(zlim[1], zlim[2], length = n.col)
    image(1, z.label, matrix(z.label, nrow = 1), zlim, c(1, 1), zlim,
          col = var.col, axes = FALSE, xlab = "", ylab = "")
    axis(4)
    box()
  dev.off()
}
