#' Data Generators
#' 
#' Collection of longitudinal responses with possibly varying lengths and 
#' varying numbers of observations. Support is [start . . . 0 . . . end] so that
#' all trajectories are aligned at 0. Usually, zero is the intervention date.
#' 
#' @section Details:
#' Generate longitudinal data for a response variable. Trajectories start
#' at time uniformly distributed in first and end at time uniformly 
#' distributed in last. Number of observations in a trajectory is 
#' Poisson(lambda_obs). The result is a number of trajectories, all starting at 
#' time 0, with different time spans, and with independently different numbers 
#' of observations within the time spans. Each trajectory follows a randomly 
#' selected response function with added N(mean, sd) error.
#' 
#' @param PL Data structure giving generation parameters
#' @param reference A nominal response reference (for example, blood pressure
#' is near 100, which is the default)
#' @param noise Vector of length 2 giving the mean and standard deviation, 
#' c(mean, sd), of added N(mean, sd) noise.
#' @param k Number of response types (TODO works only for k = 2 and 3)
#' 
#' @section Value:
#' A data frame with one response per row and three columns: id (an integer in 
#' 1:(3*n_id), time (an iteger in 0 to last observation time), and response.
#' 
#' @export
gen_traj_data = function(PL, reference = 100, noise = c(0, abs(reference/20)),
                         k = 3)
{
  ## TODO convert to PL structure
  n_id = PL$gen_par$n_id
  lambda_obs = PL$gen_par$m_obs
  first = PL$gen_par$s_range
  last = PL$gen_par$e_range
  plots = PL$gen_par$plots
  
  ## support
  ## start . . . . . 0 . . . . end
  gendata = function(n_id, lambda_obs, first, last, k) {
    id_scale = 3 # id non-consecutive scale
    
    start = floor(runif(n_id, min = first[1], max = first[2])) # 1st observation time
    end = floor(runif(n_id, min = last[1], max = last[2])) # last observation time
    n_obs = 3 + rpois(n_id, lambda = lambda_obs) # number of observations
    type = sample(k, n_id, replace = TRUE) # curve type
    id_vec = sample.int(id_scale*n_id, n_id, replace = FALSE) # make id non-consecutive
    lapply(1:n_id,   # construct by id bunches of responses
           function(i) cbind(rep(id_vec[i], n_obs[i]),
                             response(n_obs[i], type[i], start[i], end[i], reference)))
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
                  reference*sin(pi/4 + pi*times/last[2]), # part of sin curve
                  reference*(1.5 - 1/(1 + exp(-times/last[2]*5))) # 1 - sigmoid
                  )
    cbind(times, resp + rnorm(n_obs, mean = noise[1], sd = noise[2])) # add Gaussian noise
  }

  id_list = gendata(n_id, lambda_obs, first, last, k)
  id_df = as.data.frame(do.call(rbind, id_list))
  names(id_df) = c("id", "time", "response")
  
  id_df
}

## NOTES
## There are data generation packages: most notably - simstudy includes 
## longitudinal data.
## 
