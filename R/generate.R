#' Generates data for one `id`
#' 
#' @param id
#' A unique integer.
#' @param n_obs
#' An integer number of observations to produce.
#' @param type
#' Response type, 1 is constant, 2 is a sin curve portion, and 3 is a sigmoid
#' portion.
#' @param start
#' Negative integer giving time of first observation.
#' @param end
#' Positive integer giving time of last observation.
#' @param smin
#' The smallest possible `start` value among all `id`s. Currently not used.
#' @param emax
#' The largest possible `end` value among all `id`s. Used to scale sin and
#' sigmoid support.
#' 
#' @return 
#' An `n_obs` by `3 + length(vars)` matrix with columns `id`, `true_group`,
#' `time`, and responses `vars`.
#' 
oneid = function(vars, clusters, id, n_obs, type, start, end, smin, emax) {
  id = rep(id, n_obs)
  true_group = rep(type, n_obs)
  time = c(start, sort(floor(runif(n_obs - 3, min = start, max = end))), end)
  time = c(time[time <= 0], 0, time[time > 0]) # insert 0
  
  dat=as.data.table(cbind(id,true_group,time))
  line=clusters[type,2]
  #datlist=list(clus1,clus2,clus3,clus4,clus5,clus6,clus7,clus8,clus9)
  dat=addColumns(dtDefs = datlist[[line]],dat)
  
  dat
  
}

#' gendata
#' 
#' Generates data for up to three trajectory clusters
#' 
#' @param vars
#' See parameters of {\code{\link{gen_traj_data}}}.
#' @param clusters
#' See parameters of {\code{\link{gen_traj_data}}}.
#' @param m_obs
#' See parameters of {\code{\link{gen_traj_data}}}.
#' @param s_range
#' See parameters of {\code{\link{gen_traj_data}}}.
#' @param e_range
#' See parameters of {\code{\link{gen_traj_data}}}.
#' @param min_obs
#' See parameters of {\code{\link{gen_traj_data}}}.
#' @param reference
#' See parameters of {\code{\link{gen_traj_data}}}.
#' @param noise
#' See parameters of {\code{\link{gen_traj_data}}}.
#' 
#' @details 
#' Time support of each `id` is at least `s . . . 0 . . . . e`, where `s` is in
#' `s_range` and `e` is in `e_range`. 
#' @return 
#' A list of length `sum(n_id)`, where each element is a matrix output by {\code{\link{oneid}}}.
#' 
gendata = function(vars, clusters, m_obs, s_range, e_range, min_obs) {
  n_id = clusters[, "n_id"]
  t_id = sum(n_id)
  idr = c(1, 3*t_id) # id range to pick (so non-consecutive)
  
  start = round(runif(t_id, min = s_range[1], max = s_range[2]))
  end = floor(runif(t_id, min = e_range[1], max = e_range[2]))
  n_obs = min_obs + rpois(t_id, lambda = m_obs)
  type = sample(1:length(n_id), t_id, replace = TRUE, prob = n_id/t_id)
  id_vec = sample.int(idr[2], t_id, replace = FALSE) # non-consecutive
  lapply(1:t_id,   # construct by id
         function(i) oneid(vars, clusters, id_vec[i], n_obs[i], type[i],
                           start[i], end[i], smin = min(start), emax = max(end)))
}

#' Data Generators
#' 
#' Generates a collection of longitudinal responses with possibly varying
#' lengths and varying numbers of observations. Support is
#' `start` . . . 0 . . . `end`, where
#' `start`~uniform(s_range) and `end`~uniform(e_range), so that
#' all trajectories are aligned at 0 but can start and end at different times.
#' Zero is the intervention time.
#' 
#' @section Details:
#' Generate longitudinal data for a response variable. Trajectories start
#' at time uniformly distributed in s_range and end at time uniformly 
#' distributed in e_range. Number of observations in a trajectory is 
#' Poisson(m_obs). The result is a number of trajectories, all starting at 
#' time 0, with different time spans, and with independently different numbers 
#' of observations within the time spans. Each trajectory follows a randomly 
#' selected response function with added N(mean, sd) error.
#' 
#' @param vars
#' A character vecor of variable names
#' @param clusters
#' A matrix or a data frame with `length(vars)` rows and `2 + 2*length(vars)` 
#' columns giving cluster size, noise coefficient of variation (defaults to 1),
#' variable means, and variable curve types (1 to 3), in that order.
#' @param m_obs 
#' Mean number of observation per id. Provides \code{lambda} parameter in
#' \code{\link[stats]{rpois}}.
#' @param s_range 
#' A vector of length 2, giving the min and max limits of uniformly generated
#' start observation time.
#' @param e_range 
#' A vector of length 2, giving the min and max limits of uniformly generated
#' end observation time.
#' @param min_obs 
#' Minimum number of observations in addition to zero time observation.
#' @param cv
#' Coefficient of variation on generated data. If `clusters` is supplied as a
#' matrix of data.frame, it will override this parameter.
#' 
#' @return A data table with one response per row and `3 + length(vars)`
#' columns: `id`, `true_group`, `time`, and a column for each of vars.
#'   
#' @examples
#' data = gen_traj_data(n_id = c(50, 100), m_obs = 20, s_range = c(-365, -14),
#'               e_range = c(0.5*365, 2*365))
#' head(data)
#' tail(data)
#' 
#' @importFrom stats dist rnorm rpois runif
#' @export
gen_traj_data = function(vars, clusters, m_obs, curvlist, s_range, e_range, min_obs = 3,
                         verbose = FALSE){

  if(is.numeric(clusters)) {
    if(length(clusters) == 1) {
      n_id = 2^(1:clusters)*500
    } else { n_id = clusters }
    d = length(vars)*length(n_id)
    
    curv =matrix(curvlist,nrow=length(n_id),ncol=1,byrow=FALSE)
    
    clusters = cbind(n_id, curv)
    colnames(clusters) = c("n_id", paste0("curve"))
  } else {
    if(class(clusters) == "matrix" | class(clusters == "data.frame")) {
      if(ncol(clusters) != 2 + 2*nrow(clusters)) 
        stop("gen_traj_data: clusters parameter has invalid dimensions")
    } else {
      stop("gen_traj_data: clusters parameter invalid class")
    }
  }
  if(verbose) print(clusters)
  
  id_list = gendata(vars, clusters, m_obs, s_range, e_range, min_obs)
  id_mat = do.call(rbind, id_list)
  Datalist=data.table::data.table(id_mat)
  #}
 # return(Datalist)
}
                                              
