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
#' @param n_id Number of trajectories
#' @param lambda_obs Poisson mean for number of observations in a trajectory
#' over 3. That is, the mean is lambda_obs + 3.
#' @param first Vector of length 2 giving first observation time range, 
#' c(min, max). First observation is uniformly distributed over this range. 
#' @param last Vector of length 2 giving last observation time range, 
#' c(min, max). Last observation is uniformly distributed over this range. 
#' @param reference A nominal response reference (for example, blood pressure
#' is near 100, which is the default)
#' @param noise Vector of length 2 giving the mean and standard deviation, 
#' c(mean, sd), of added N(mean, sd) noise.
#' @param plots Number of plots of randomly selected trajectories to check
#' generated data
#' 
#' @section Value:
#' A data frame with one response per row and three columns: id (an integer in 
#' 1:n_id), time (an iteger in 0 to last observation time), and response.
#' 
#' @export
gen_traj_data = function(n_id, lambda_obs, first = c(-50, -10), last = c(50,100),
                         reference = 100,
                         noise = c(0, abs(reference/20)), k = 3, plots = 0)
{
  # n_id = PL$gen_par$n_id
  # lambda_obs = PL$gen_par$lambda_obs
  # first = PL$gen_par$first
  # last = PL$gen_par$last
  # plots = PL$gen_par$plots
  
  ## support
  ## start . . . . . 0 . . . . end
  gendata = function(n_id, lambda_obs, first, last, k) {
    start = floor(runif(n_id, min = first[1], max = first[2])) # 1st observation time
    end = floor(runif(n_id, min = last[1], max = last[2])) # last observation time
    n_obs = 3 + rpois(n_id, lambda = lambda_obs) # number of observations
    type = sample(k, n_id, replace = TRUE) # curve type
    id_vec = sample.int(3*n_id, n_id, replace = FALSE) # make id non-consecutive
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
                  reference*sin(pi/4 + 0.5*pi*times/last[2]), # part of sin curve
                  reference*(1.5 - 1/(1 + exp(-times/last[2]*5))) # 1 - sigmoid
                  )
    cbind(times, resp + rnorm(n_obs, mean = noise[1], sd = noise[2])) # add Gaussian noise
  }

  id_list = gendata(n_id, lambda_obs, first, last, k)
  id_df = as.data.frame(do.call(rbind, id_list))
  names(id_df) = c("id", "time", "response")

  ## plot to check result
  if(plots > 0) {
    require(ggplot2)
    iplot = sample(unique(id_df$id), plots)
    print(ggplot(dplyr::filter(id_df, id %in% iplot), aes(x = time, y = response)) +
            facet_wrap(~ id) + geom_point())
  }
  
  id_df
}

## NOTES
## There are data generation packages: most notably - simstudy includes 
## longitudinal data.
## 
