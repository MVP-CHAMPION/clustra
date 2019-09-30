#' Data Generators
#' 
#' Collection of longitudinal responses with possibly varying lengths and 
#' varying numbers of observations.
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
#' @param n_id Number of trajectories
#' @param m_obs Poisson mean for number of observations in a trajectory
#' @param s_range Vector of length 2 giving last observation time range, 
#' c(min, max). Lengths are uniformly distributed over this range. 
#' @param e_range Vector of length 2 giving last observation time range, 
#' c(min, max). Lengths are uniformly distributed over this range. 
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
gen_long_data = function(n_id, m_obs, s_range = c(0, 0), 
                         e_range = c(50,100), center = 100, 
                         noise = c(0, center/20), n_typ = 3, plots = 0)
{
  ## support
  gendata = function(n_id, m_obs, s_range, e_range, n_typ) {
    start = runif(n_id, min = s_range[1], max = s_range[2]) # 1st observation time
    end = runif(n_id, min = e_range[1], max = e_range[2]) # last observation time
    n_obs = 1 + rpois(n_id, lambda = m_obs) # number of observations
    type = sample(n_typ, n_id, replace = TRUE) # curve type
    id_vec = sample.int(3*n_id, n_id, replace = FALSE) # make id non-consecutive
    lapply(1:n_id,   # construct by id bunches of responses
           function(i) cbind(rep(id_vec[i], n_obs[i]),
                             response(n_obs[i], type[i], start[i], end[i], center)))
  }
  ## response
  response = function(n_obs, type, start, end, center) {
    times = c(start, sort(runif(n_obs - 1))*end)
    ## TODO add user supplied functions
    resp = switch(type, # 1 = constant, 2 = sin, 3 = sigmoid
                  rep(center, length(times)), # constant
                  center*sin(pi/4 + 0.5*pi*times/e_range[2]), # part of sin curve
                  center*(1.5 - 1/(1 + exp(-times/e_range[2]*5))) # 1 - sigmoid
                  )
    cbind(times, resp + rnorm(n_obs, mean = noise[1], sd = noise[2])) # add Gaussian noise
  }

  id_list = gendata(n_id, m_obs, s_range, e_range, n_typ)
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
