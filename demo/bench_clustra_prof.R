## profile clustra
set.seed(123473)
data = gen_traj_data(200000, m_obs = 45, s_range = c(-50, -10), e_range = c(3*365, 10*365), reference = 100)
data$id = as.numeric(factor(data$id))
set.seed(123473)
Rprof()
cl = clustra(data, k = 3, fp = list(maxdf = 30, iter = 10, starts = 5, idperstart = 20, retry_max = 3), cores = c(e_mc = 1, m_mc = 1, nthreads = 1, blas = 1), verbose = TRUE)
Rprof(NULL)
summaryRprof()

library(profvis)
library(bench)
set.seed(123477)
profvis({
  group = sample(3, length(unique(data$id)), replace = TRUE)
  fp = list(maxdf = 30, iter = 1, starts = 5, idperstart = 20, retry_max = 3)
  data$group = group[data$id] # expand group to all data
  f = trajectories(data, 3, group, fp = fp, cores = c(e_mc = 1, m_mc = 1, nthreads = 1, blas = 1), verbose = TRUE)
  if(!is.null( (er = xit_report(f, fp)) )) print(er)
})

