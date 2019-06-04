## generate longitudinal blood pressure data
#' @param n_id Number of trajectories
#' @param m_obs Poisson mean for number of observations per trajectory 
#' @export
gen_bp_data = function(n_id, m_obs)
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

  gendata = function(n_id, m_obs, n_typ = 3, t_max = 365*10) {
    id_list = vector("list", n_id)
    for(i in seq(1, n_id)) {
      n_obs = 1 + rpois(1, lambda = m_obs)
      type = sample(n_typ, 1)
      ti_max = runif(1, min = 600, max = t_max)
      id_list[[i]] = cbind(rep(i, n_obs), response(n_obs, type, ti_max))
    }
    id_list
  }

  id_list = gendata(n_id, m_obs)
  id_df = as.data.frame(do.call(rbind, id_list))
  names(id_df) = c("id", "age", "sbp")

  id_df

  ## plot to check result
  ## for(id in id_list) plot(c(0, id[, 2]), c(0, id[, 3]))
}