#' Data Generators
#' 
#' Generates a collection of longitudinal responses with possibly varying
#' lengths and varying numbers of observations. Support is
#' [start . . . 0 . . . end], where
#' start~uniform(s_range) and end~uniform(e_range), so that
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
#' @param n_id 
#' Number of id to generate.
#' @param m_obs 
#' Mean number of observation per id. Provides \code{lambda} parameter in
#' \code{\link[stats]{rpois}}.
#' @param s_range 
#' A vector of length 2, giving the min and max limits of uniformly generated
#' start observation time.
#' @param e_range 
#' A vector of length 2, giving the min and max limits of uniformly generated
#' end observation time.
#' @param reference 
#' A nominal response value (for example, blood pressure is near 100, which
#' is the default)
#' @param noise 
#' Vector of length 2 giving the *mean* and *sd* of added N(mean, sd) noise.
#' @param k 
#' Number of response types (So far works only for k = 2 and 3)
#' @param min_obs 
#' Minimum number of observations in addition to zero time observation.
#' 
#' @return A data frame with one response per row and three columns:
#'   \code{id}, \code{time}, \code{response}, and \code{group}.
#' @examples
#' set.seed(123)
#' data = gen_traj_data(n_id = 20, m_obs = 10, s_range = c(-50, -10),
#'               e_range = c(2*365, 5*365))
#' head(data)
#' 
#' @importFrom stats dist rnorm rpois runif
#' @export
gen_traj_data = function(n_id, m_obs, s_range, e_range, reference = 100,
                         noise = c(0, abs(reference/20)), k = 3, min_obs = 3)
{
  ## support
  ## s_range . . . 0 . . . . e_range
  gendata = function(n_id, m_obs, s_range, e_range, k) {
    idr = c(1, 3*n_id) # id range to pick (so non-consecutive)
    
    start = round(runif(n_id, min = s_range[1], max = s_range[2])) # start time
    end = floor(runif(n_id, min = e_range[1], max = e_range[2])) # end time
    n_obs = min_obs + rpois(n_id, lambda = m_obs) # number of observations
    type = sample(k, n_id, replace = TRUE) # curve type
    id_vec = sample.int(idr[2], n_id, replace = FALSE) # non-consecutive
    lapply(1:n_id,   # construct by id
           function(i) oneid(id_vec[i], n_obs[i], type[i], start[i], end[i],
                                reference))
  }
  ## response
  oneid = function(id, n_obs, type, start, end, reference) {
    id = rep(id, n_obs)
    true_group = rep(type, n_obs)
    time = c(start,
              start + sort(floor(runif(n_obs - 3, min = start, max = end))),
              end)
    time = c(time[time <= 0], 0, time[time > 0]) # insert 0

    ## TODO add user supplied functions
    response = switch(type, # 1 = constant, 2 = sin, 3 = sigmoid
                  rep(reference, length(time)), # constant
                  reference*sin(pi/4 + pi*time/e_range[2]), # part of sin curve
                  reference*(1.5 - 1/(1 + exp(-time/e_range[2]*5))) # 1 - sigmoid
                  ) + rnorm(n_obs, mean = noise[1], sd = noise[2])
    cbind(id, time, response, true_group)
  }

  id_list = gendata(n_id, m_obs, s_range, e_range, k)
  id_mat = do.call(rbind, id_list)
  id_dt = data.table::data.table(id_mat)

  id_dt
}
