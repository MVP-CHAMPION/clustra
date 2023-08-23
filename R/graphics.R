#' Timing function
#'
#' @param ltime
#' Result of last call to deltime. 
#' @param text 
#' Text to display along with elapsed time since \code{ltime}.
#' @param units
#' Logical. If TRUE, print units
#' @param nl
#' Logical. If TRUE, a newline is added at the end.
#'
#' @return
#' "elapsed" component of current \code{\link{proc.time}}.
#' 
#' @export
deltime = function(ltime = proc.time()["elapsed"], text = NULL, units = FALSE,
                   nl = FALSE) {
  time = proc.time()["elapsed"]
  if(!is.null(text)) {
    x = round(difftime(time, ltime), 2)
    if(units) units = paste0(" ", attr(x, "units"))
    else units = ""
    if(nl) nl = "\n"
    else nl = ""
    cat(paste0(text, x, units, nl))
  }
  invisible(time)
}

#' Plots a sample of ids in a small mutiples layout
#' 
#' @param dat 
#' A data frame with a few id trajectories to plot.
#' @param layout 
#' The small multiples layout as c(rows, columns).
#' @param sample
#' If zero, all data in dat are displayed. If >0 a sample of that many data 
#' points from dat are displayed.
#' @param group
#' If not NULL, a character string giving the variable name in data that should
#' color the data points.
#' 
#' @return 
#' Invisibly returns the number of trajectories plotted.
#' 
#' @importFrom graphics split.screen screen close.screen
#' @export
plot_sample = function(dat, layout = c(3,3), sample = prod(layout),
                       group = NULL) {
  id = NULL # for data.table R CMD check
  
  ids = unique(dat$id)
  if(sample) {
    ids = sample(ids, sample)
    dat = dat[id %in% ids]
  }
  xrange = range(dat$time)
  yrange = range(dat$response)
  xat = pretty(xrange)
  yat = pretty(yrange)
  
  m = matrix(c(0.1, .95, 0.1, .95), nrow = 1)
  screen = 1
  for(i in 1:length(ids)) {
    idat = dat[id == ids[i]]
    if((screen %% prod(layout)) == 1) { #} && i > prod(layout)) {
      screen = 1
      close.screen(all.screens = TRUE)
      sc = split.screen(m)
      screen(1)
      sc = split.screen(layout)
    }
    screen((screen = screen + 1))
    opar = par(mar = c(0, 0, 0, 0))
    if(is.null(group))
      plot(idat$time, idat$response, axes = FALSE, type = "p", xlim = xrange,
           ylim = yrange, pch = 20)
    else
      plot(idat$time, idat$response, axes = FALSE, type = "p", xlim = xrange,
           ylim = yrange, col = idat[[group]] + 1, pch = 20)
    text(x = xrange[2], y = yrange[2], labels = sprintf("id=%s", ids[i]), 
         adj = c(1, 1))
    box()
    if((screen - 1) %% layout[2] == 1 || layout[2] == 1) axis(2, at = yat)
    if((screen - 1) > (layout[1] - 1)*layout[2] ||
       i > (length(ids) - layout[2])) axis(1, at = xat)
  }
  close.screen(all.screens = TRUE)
  invisible(length(ids))
}

#' Plots data and smooths from `clustra` output or internally from within
#' `start_groups`
#'
#' @param data
#' The data. If after `clustra` run, it includes resulting clusters as group.
#' @param fits 
#' The `cl$tps` component of `clustra` output or internal `start_groups` fits. 
#' If fits are supplied and `select.data` is NULL, the data is colored by 
#' clusters. If NULL, or if `select.data` is not NULL, the data is black points.
#' @param max.data 
#' The maximum number of data points to plot (defaults to 10,000). If zero,
#' no points are plotted (overrides select.data). Use `Inf` value to plot all
#' points.
#' @param select.data 
#' Either NULL or a list of length k, each element a data.frame (like data)
#' with time and response components. The select.data points will be
#' highlighted with cluster colors on the plot. This is used internally in
#' `start_groups` function to show the selected starting points. In this case,
#' also the fits parameter can contain TPS fits to the starting points.
#' @param group
#' Character variable name in `data` to color the clusters. A NULL will produce
#' a b&w point plot.
#' 
#' @importFrom graphics lines points
#' @export
plot_smooths = function(data, fits = NULL, max.data = 20000, 
                        select.data = NULL, group = "group", ...) {
  k = length(fits)
  xrng = range(data$time)
  ptime = seq(xrng[1], xrng[2], length.out = 100)
  yrng = range(data$response)
  sdt = data
  if(max.data && nrow(data) > max.data) # reduce to max.data
    sdt = data[sample(nrow(data), max.data)]
  
  ## data plotting
  if(max.data) {  # a zero will skip data plotting
    if(is.null(group)) {
      plot(sdt$time, sdt$response, pch = ".",
           xlab = "time", ylab = "response", ...)
      } else {
        plot(sdt$time, sdt$response, pch = ".", 
             col = as.numeric(sdt[[group]]) + 1,
             xlab = "time", ylab = "response", ...)
      }
    if(!is.null(select.data)) { # big dot color for selected data
      for(i in 1:length(select.data))
        points(select.data[[i]]$time, select.data[[i]]$response, pch = 20, 
               col = i + 1)
    }
  } else { # just axes and no data plot
    plot(x = sdt$time, y = sdt$response, type = "n", ...)
  }

  ## tps plotting  
  if(!is.null(fits)) { # plot color fitted lines
    for(i in 1:k) {
      pred = predict(fits[[i]], newdata = data.frame(time = ptime))
      lines(ptime, pred, col = i + 1, lwd = 2)
    }
  }
}

#' Plots a list item, a silhouette, from the result of `clustra_sil` along 
#' with the average silhouette value. Typically used via 
#' `lapply(list, plot_silhouette)`
#' @param sil
#' A data frame that is a list item returned by `clustra_sil`.
#' 
#' @return
#' Returns invisibly the average silhouette value.
#' 
#' @importFrom graphics barplot legend abline text
#' @export
plot_silhouette = function(sil) {
  k = length(levels(sil$cluster))
  barplot(height = sil$silhouette, col = as.numeric(sil$cluster) + 1,
          border = as.numeric(sil$cluster) + 1)
  legend("topright", legend = 1:k, fill = 1:k + 1)
  pct = round(mean(sil$silhouette), 2)
  abline(h = pct, lty = 5, col = "red")
  mtext(paste("Average Width", pct), col = "red")
#  text(x = nrow(sil), y = pct, labels = sprintf("Average Width %s", pct), 
#       adj = c(1, 1))
  invisible(pct)
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
