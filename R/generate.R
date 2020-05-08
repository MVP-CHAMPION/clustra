#' Data Generators
#' 
#' Collection of longitudinal responses with possibly varying lengths and 
#' varying numbers of observations. Support is [st . . . 0 . . . en] so that
#' all trajectories are aligned at 0. Zero is the intervention time.
#' 
#' @section Details:
#' Generate longitudinal data for a response variable. Trajectories start
#' at time uniformly distributed in st and end at time uniformly 
#' distributed in en. Number of observations in a trajectory is 
#' Poisson(lambda_obs). The result is a number of trajectories, all starting at 
#' time 0, with different time spans, and with independently different numbers 
#' of observations within the time spans. Each trajectory follows a randomly 
#' selected response function with added N(mean, sd) error.
#' 
#' @param n_id See .clustra_env
#' @param lambda_obs See .clustra_env
#' @param st See .clustra_env
#' @param en See .clustra_env
#' @param plots See .clustra_env
#' @param reference A nominal response reference (for example, blood pressure
#'   is near 100, which is the default)
#' @param noise Vector of length 2 giving the mean and sd of added N(mean, sd)
#'   noise.
#' @param k Number of response types (TODO works only for k = 2 and 3)
#' @param min_obs Minimum number of observations in addition to zero time
#'   observation.
#' 
#' @return A data frame with one response per row and three columns:
#'   `id` (an integer in 1:(3*n_id), `time` (an iteger in `st[1]` to `en[2]`
#'   observation time), and `response`.
#' 
#' @export
gen_traj_data = function(n_id = clustra_env("gen$n_id"),
                         lambda_obs = clustra_env("gen$m_obs"),
                         st = clustra_env("gen$s_range"),
                         en = clustra_env("gen$e_range"),
                         plots = clustra_env("gen$plots"),
                         reference = 100, noise = c(0, abs(reference/20)),
                         k = 3, min_obs = 3)
{
  ## support
  ## st . . . 0 . . . . en
  gendata = function(n_id, lambda_obs, st, en, k) {
    idr = c(1, 3*n_id) # id range to pick (so non-consecutive)
    
    start = round(runif(n_id, min = st[1], max = st[2])) # start observe time
    end = floor(runif(n_id, min = en[1], max = en[2])) # end observe time
    n_obs = min_obs + rpois(n_id, lambda = lambda_obs) # number of observations
    type = sample(k, n_id, replace = TRUE) # curve type
    id_vec = sample.int(idr[2], n_id, replace = FALSE) # non-consecutive
    lapply(1:n_id,   # construct by id bunches of responses
           function(i) cbind(rep(id_vec[i], n_obs[i]),
                             response(n_obs[i], type[i], start[i], end[i], 
                                      reference)))
  }
  ## response
  response = function(n_obs, type, start, end, reference) {
    times = c(start,
              start + sort(floor(runif(n_obs - 3, min = start, max = end))),
              end)
    times = c(times[times <= 0], 0, times[times > 0]) # insert 0
    ## TODO add user supplied functions
    resp = switch(type, # 1 = constant, 2 = sin, 3 = sigmoid
                  rep(reference, length(times)), # constant
                  reference*sin(pi/4 + pi*times/en[2]), # part of sin curve
                  reference*(1.5 - 1/(1 + exp(-times/en[2]*5))) # 1 - sigmoid
                  )
    cbind(times, resp + rnorm(n_obs, mean = noise[1], sd = noise[2])) # add Gaussian noise
  }

  id_list = gendata(n_id, lambda_obs, st, en, k)
  id_df = as.data.frame(do.call(rbind, id_list))
  names(id_df) = c("id", "time", "response")
  
  id_df
}
