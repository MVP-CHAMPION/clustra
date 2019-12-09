## This is the R translation of trajectories v2.sas code
# library(data.table) # possibly for reading data
library(gratia) # from GitHub: "gavinsimpson/gratia"
library(parallel)
library(mgcv)
library(proftools)
library(openblasctl)

source("R/deltime.R")
a0 = a = deltime()

## Function to fit thin plate spline (tps) to a group with gam from the mgcv
## package. Crossvalidation does not know about zero weights, resulting in 
## different smoothing parameters, so subset is used and a separate prediction
## function.
## TODO Add best of multiple random starts based on minimum deviance.
## TODO Add uncertainty moderated outlier detection. Do we need to use credible
##      regions or is se.fit enough??
## TODO Explore relationship with bam (a "large" version of gam). Initial tries
##      took loner to compute.
## TODO Add Rand index determination of k, explore AIC and BIC. Lots of parallel
## opportunities.
## TODO Configure with gam-bam choice for benchmaarking #####################
tps_g = function(g, dat, maxdf) {
  tps = mgcv::bam(response ~ s(time, k = maxdf), data = dat, 
                  subset = dat$group == g, discrete = TRUE)
#  tps = mgcv::gam(response ~ s(time, k = maxdf), data = dat, 
#                  subset = dat$group == g)
  if(class(tps)[1] == "try-error") print(attr(tps, "condition")$message)
  pred = predict(tps, newdata = dat, type = "response", se.fit = TRUE)
  pred$tps = tps
  pred
}

## Loss functions for each id to a group tps center
mse_g = function(g, tps, dat) # mean squared error
  as.numeric(tapply((dat$response - tps[[g]]$fit)^2, dat$id, mean))
mxe_g = function(g, tps, dat) # maximum error
  as.numeric(tapply(abs(dat$response - tps[[g]]$fit), dat$id, max))
## Fn to compute points outside CI
## Mean distance in probability (a multi-Mahalanobis distance)
mpd_g = function(g, tps, dat) {
  ## TODO 
}


#' Cluster longitudinal trajectories of a response variable. Trajectory means
#' are splines fit to all ids in a cluster.
#' @param dat Data frame with response measurements, one per observation. Column
#' names are id, time, response.
#' @param ng Number of clusters
#' @param iter Maximum number of iterations.
#' @param maxdf Maximum degrees of freedom for trajectory spline smooths.
#' @param plot Plot clustered data with superimposed trajectory spline smooths.
#' @export
trajectories = function(dat, ng, iter = 10, maxdf = 50, plot = FALSE) {
  ## get number of unique id
  if(plot) require(ggplot2)
  n_id = length(unique(dat$id)) 
  
  ## start with random group assignments
  group = sample(ng, n_id, replace = TRUE)
  dat$id = as.numeric(factor(dat$id)) # removes assumption id are sequential
  dat$group = group[dat$id] # expand group assignment to all responses
    
  ## EM algorithm to cluster ids into ng groups
  ## iterate fitting a thin plate spline center to each group (M-step)
  ##         regroup each id to nearest tps center (E-step)
  a = deltime(a, "Starting iteration")
  for(i in 1:iter) {
    ## M-step:
    ##   fit tp spline centers for each group separately
    pd = profileExpr({
      tps = lapply(1:ng, tps_g, dat = dat, maxdf = maxdf)
    })
    #      tps = mclapply(1:ng, tps_g, dat = dat, maxdf = maxdf, mc.cores = 1)
      a = deltime(a, "M-step")

    ## E-step:
    ##   compute loss of each id to each tp spline
    loss = mclapply(1:ng, mse_g, tps = tps, dat = dat, mc.cores = 1)
    ## get mse-nearest tp spline for each id
    new_group = apply(do.call(cbind, loss), 1, which.min)
    a = deltime(a, "E-step")
    
    ## Check if done or degenrate
    changes = sum(new_group != group)
    counts = tabulate(new_group)
    deviance = sum(unlist(lapply(1:ng, function(g) deviance(tps[[g]]$tps))))
    cat("Iteration:", i, "Changes:", changes, "Counts:", counts,
        "Deviance:", deviance, "\n")
    if(length(counts) < ng | any(counts == 0)) {
      cat("Empty group encountered. Restarting on next iteration.\n")
      new_group = sample(ng, n_id, replace = TRUE)
      counts = tabulate(new_group)
    }
    group = new_group
    dat$group = as.factor(group[dat$id])
    if(changes == 0) break
    
    ## add outlier treatment as probability of belonging?
    ## use predict.gam(), probably for each group separately? weights??
  }
  a = deltime(a, "Done iteration")
  
  if(plot) print(
    ggplot(dat, aes(time, response, color = group, group = group)) +
      geom_point(alpha = 0.1) +
      stat_smooth(method = "gam", formula = y ~ s(x, k = maxdf), size = 1) +
      theme_bw())
  a = deltime(a, "Done plots")
  print(funSummary(pd))
  print(funSummary(pd, srclines = FALSE))
  print(callSummary(pd))
  print(srcSummary(pd))
  print(hotPaths(pd, total.pct = 10.0))
  plotProfileCallGraph(pd)
  plotProfileCallGraph(pd, focus = "mgcv::bam")
  flameGraph(pd)
  calleeTreeMap(pd)
  
  list(group = group, dat_group = dat$group, tps = tps)
}

sessionInfo()
source("R/generate.R")
set.seed(90)
dat = gen_long_data(n_id = 2000000, m_obs = 25, e_range = c(365*3, 365*10),
                    plots = FALSE)
a = deltime(a, paste0("Data (", paste(dim(dat), collapse = ","), ") generated"))

openblas_set_num_threads(1)
  f = trajectories(dat = dat, ng = 3, iter = 1, maxdf = 50, plot = FALSE)
## pd_gam = filterProfileData(pd, select = "gam.setup")
## funSummary(pd_gam)
## funSummary(pd_gam, srclines = FALSE)
## callSummary(pd_gam)
## srcSummary(pd_gam)
## hotPaths(pd_gam, total.pct = 10.0)
## plotProfileCallGraph(pd_gam)
## flameGraph(pd_gam)
## calleeTreeMap(pd_gam)

dat$group = f$dat_group

# source("R/benchmark.R")
# set.seed(90)
# dat = gen_long_data(n_id = 1000, m_obs = 25, e_range = c(365*3, 365*10),
#                     plots = 20)
# bench_cores(FUN = trajectories, dat = dat, ng = 3, iter = 20, maxdf = 50,
#             plot = FALSE, max2 = 1, reps = 1)

a = deltime(a0, "Total time")

