#' Generates data for one `id`
#' 
#' @param id
#' A unique integer.
#' @param n_obs
#' An integer number of observations to produce.
#' @param type
#' Response type, 1 is constant, 2 is a sin curve portion, and 3 is a sigmoid
#' portion.
#' @param intercept 
#' Used to set response at `smin` time value (not 0) and shift all responses 
#' accordingly.
#' @param start
#' Negative integer giving time of first observation.
#' @param end
#' Positive integer giving time of last observation.
#' @param smin
#' The smallest possible `start` value among all `id`s. Used to align with 
#' intercept and then dropped.
#' @param emax
#' The largest possible `end` value among all `id`s. Used to scale sin and
#' sigmoid support.
#' @param noise
#' Standard deviation of zero mean Gaussian noise added to response functions.
#' @return 
#' An `n_obs` by 4 matrix with columns `id`, `time`, `response`, `true_group`.
#' 
oneid = function(id, n_obs, type, intercept, start, end, smin, emax, noise) {
  id = rep(id, n_obs)
  true_group = rep(type, n_obs)
  time = c(smin, start, sort(floor(runif(n_obs - 3, min = start, max = end))), end)
  time = c(time[time <= 0], 0, time[time > 0]) # insert 0
  
  base = switch(type, # 1 = constant, 2 = sin, 3 = sigmoid
                    100*rep(0, length(time)), # constant
                    100*sin(2*pi/3 + pi*time/emax), # part of sin curve
                    100*(1 - 1/(1 + exp(-time/emax*5))) # 1 - sigmoid
  )
  base = intercept + (base - base[1])
  response = base[-1] + rnorm(n_obs, mean = noise[1], sd = noise[2])
  time = time[-1]
  cbind(id, time, response, true_group) # removed smin time used to align
}

#' gendata
#' 
#' Generates data for up to three trajectory clusters
#' 
#' @param n_id
#' See parameters of {\code{\link{gen_traj_data}}}.
#' @param m_obs
#' See parameters of {\code{\link{gen_traj_data}}}.
#' @param types
#' See parameters of {\code{\link{gen_traj_data}}}.
#' @param intercepts
#' See parameters of {\code{\link{gen_traj_data}}}.
#' @param s_range
#' See parameters of {\code{\link{gen_traj_data}}}.
#' @param e_range
#' See parameters of {\code{\link{gen_traj_data}}}.
#' @param min_obs
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
gendata = function(n_id, types, intercepts, m_obs, s_range, e_range, min_obs,
                   noise) {
  idr = c(1, 3*sum(n_id)) # id range to pick (so non-consecutive)
  t_id = sum(n_id)
  
  start = round(runif(t_id, min = s_range[1], max = s_range[2]))
  end = floor(runif(t_id, min = e_range[1], max = e_range[2]))
  n_obs = min_obs + rpois(t_id, lambda = m_obs)
  type = sample(1:length(types), t_id, replace = TRUE, prob = n_id/t_id)
  id_vec = sample.int(idr[2], t_id, replace = FALSE) # non-consecutive
  rlist = lapply(1:t_id,   # construct by id
                 function(i) oneid(id_vec[i], n_obs[i], types[type[i]], 
                                   intercepts[type[i]], start[i], end[i], 
                                   smin = min(start), emax = max(end), noise))
  return(rlist)
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
#' of observations within the time spans. Each trajectory follows one of three 
#' possible response functions possibly with a different mean and with added
#' N(mean, sd) error.
#' 
#' @param n_id 
#' Vector whose length is the number of clusters, giving the number of id's to
#' generate in each cluster.
#' @param types 
#' A vector of integers from `c(1, 2, 3)` of same length as n_id, indicating 
#' curve type: constant, sine portion, sigmoid portion, respectively.
#' @param intercepts
#' A vector of first responses at minimum time for the curve base vectors of 
#' same length as n_id. Each `type`-`intercept` combination should be unique for 
#' unique clusters.
#' @param m_obs 
#' Mean number of observation per id. Provides \code{lambda} parameter in
#' \code{\link[stats]{rpois}}.
#' @param s_range 
#' A vector of length 2, giving the min and max limits of uniformly generated
#' start observation time.
#' @param e_range 
#' A vector of length 2, giving the min and max limits of uniformly generated
#' end observation time.
#' @param noise 
#' Vector of length 2 giving the *mean* and *sd* of added N(mean, sd) noise.
#' @param min_obs 
#' Minimum number of observations in addition to zero time observation.
#' 
#' @return A data table with one response per row and four columns:
#'   `id`, `time`, `response`, and `true_group`.
#'   
#' @examples
#' data = gen_traj_data(n_id = c(50, 100), types = c(1, 2), 
#'   intercepts = c(100, 80), m_obs = 20, s_range = c(-365, -14), 
#'   e_range = c(0.5*365, 2*365))
#' head(data)
#' tail(data)
#' 
#' @importFrom stats dist rnorm rpois runif
#' @export
gen_traj_data = function(n_id, types, intercepts, m_obs, s_range, e_range,
                         noise = c(0, abs(mean(intercepts)/20)), min_obs = 3)
{
  id_list = gendata(n_id, types, intercepts, m_obs, s_range, e_range, min_obs, 
                    noise)
  id_mat = do.call(rbind, id_list)
  data = data.table::data.table(id_mat)
  return(data)
}
