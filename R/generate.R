## generate longitudinal blood pressure data
#' @param n_id Number of trajectories
#' @param m_obs Poisson mean for number of observations in a trajectory
#' @param lrange Vector c(min, max) range for last observation time generation 
#' @export
gen_bp_data = function(n_id, m_obs, lrange = c(365*3, 365*10), plots = FALSE)
{
  ## response functions
  response = function(n_obs, type = 1, ti_max,
                      center = 100) {
    times = c(0, sort(runif(n_obs - 1)))*ti_max
    resp = switch(type,
                  rep(center, length(times)), # constant
                  center*sin(pi/4 + pi*times/ti_max/2), # part of sin curve
                  tanh((times - ti_max/2)*2*pi/ti_max)*(center/5) + center # tanh
                  )
    cbind(times, resp + rnorm(n_obs, sd = center/20)) # add Gaussian noise
  }

  ## support
  gendata = function(n_id, m_obs, n_typ = 3, t_min = lrange[1], t_max = lrange[2]) {
    ti_max = runif(n_id, min = t_min, max = t_max) # last observation
    n_obs = 1 + rpois(n_id, lambda = m_obs) # number of observations
    type = sample(n_typ, n_id, replace = TRUE) # curve type
    lapply(1:n_id, function(i) cbind(rep(i, n_obs[i]), response(n_obs[i], type[i], ti_max[i])))
  }

  id_list = gendata(n_id, m_obs)
  id_df = as.data.frame(do.call(rbind, id_list))
  names(id_df) = c("id", "age", "sbp")

  ## plot to check result
  if(plots) {
    iplot = sample(unique(id_df$id), plots)
    print(ggplot(filter(id_df, id %in% iplot), aes(x = age, y = sbp)) +
            facet_wrap(~ id) + geom_point())
  }
  
  id_df
}