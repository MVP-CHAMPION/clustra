## Various functions to evaluate results of clustering, like RandIndex
## computation and plots of spline fits.
##
#' Plots spline fit from bam, possibly with data
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

## Use RandIndex to evaluate number of clusters
randplot = function(results) {
  library(MixSim)
  nr = length(results)
  rand_mat = vector("list", nr*(nr - 1)/2)
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
      rand_mat[[irm]] = data.frame(i.K = i.K, i.R = i.R, j.K = j.K, j.R = j.R, 
                                   Rand = Rand, adjRand = adjRand, Fow = Fow, Mir = Mir)
    }
  }
  rand_mat = do.call(rbind, rand_mat)
  rand_plot(rand_mat)
}

## Rand index matrix plot from Technometrics paper
## Sorts replicates within cluster K
## Author: Wei-Chen Chen
## 
#' @param randata A data frame with columns 
rand_plot = function(randata) {
  library(RColorBrewer)

  n.col = max(randata$i.K)*max(randata$i.R)
  randata <- as.data.frame(randata)
  K.vec <- unique(unlist(randata[, c("i.K", "j.K")]))
  K.max <- max(K.vec)
  K.len = length(K.vec)
  unique.R <- unique(unlist(randata[, c("i.R", "j.R")]))
  R.vec <- max(unique.R)
  R.max = max(R.vec)

  x <- 1:((K.max - 1) * R.vec)
  y <- x
  z <- array(NA, c(max(x), max(x)))
  diag(z) <- 1
  var.z <- randata$adjRand
  for(i in 1:nrow(randata)){
    z.i <- (randata$i.K[i] - 2) * R.vec + randata$i.R[i]
    z.j <- (randata$j.K[i] - 2) * R.vec + randata$j.R[i]
    z[z.i, z.j] <- z[z.j, z.i] <- var.z[i]
  }
  lim <- range(x)
  zlim <- range(c(randata$adjRand, 1))
  var.col <- colorRampPalette(brewer.pal(9, "YlOrRd"))(n.col)

## Reorder within K for vis effect
  for(i in 1:K.len) { 
#    id.z <- order(ret.ic[ret.ic[, 1] == (i + 1), 3], decreasing = TRUE)
    id.z <- order(z[(i - 1) * R.max + (1:R.max), (i - 1) * R.max + 1], decreasing = TRUE)
    tmp.z <- z[(i - 1) * R.max + (1:R.max), ]
    z[(i - 1) * R.max + (1:R.max), ] <- tmp.z[id.z,]
    tmp.z <- z[, (i - 1) * R.max + (1:R.max)]
    z[, (i - 1) * R.max + (1:R.max)] <- tmp.z[, id.z]

##    cat("K = ", i + 1, ", R order =", (1:10)[id.z], "\n")
}

  bg <- "transparent"
  fg <- "black"

  pdf("all_sort1.pdf", width = 6, height = 6, bg = bg)
    par(bg = bg, fg = fg, col = fg, col.axis = fg,
        col.lab = fg, col.main = fg, col.sub = fg)
        #    cex.main = 2.0, cex.lab = 2.0, cex.axis = 1.5)
        #layout(matrix(1:2, nrow = 1), width = c(7, 2))

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

    z.label <- seq(zlim[1], zlim[2], length = n.col)
    image(1, z.label, matrix(z.label, nrow = 1), zlim, c(1, 1), zlim,
          col = var.col, axes = FALSE, xlab = "", ylab = "")
    axis(4)
    box()
  dev.off()
}
