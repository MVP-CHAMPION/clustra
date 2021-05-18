## profile clustra
library(clustra)
set.seed(123)
data = gen_traj_data(320000, m_obs = 45, s_range = c(-50, -10), e_range = c(3*365, 10*365), reference = 100)
library(data.table)
data = as.data.table(data)
data = data[order(id)]
data[, id := as.numeric(factor(id))]
group = sample(3, data[, length(unique(id))], replace = TRUE)
fp = list(maxdf = 30, iter = 10, starts = 5, idperstart = 20, retry_max = 3)
data[,group:=group[id]]
set.seed(1234)
f = trajectories(data, 3, group, fp = fp, cores = c(e_mc = 1, m_mc = 1, nthreads = 1, blas = 1), verbose = TRUE)

## Above works fast with data.table 23 sec/iteration on 320K id's. But it is
## probably not threadsafe precluding the use of mclapply! Is data.table
## threaded for by= operations?? Multihthreading over id's would be great!

set.seed(1234)
Rprof()
cl = clustra(data, k = 3, fp = list(maxdf = 30, iter = 10, starts = 12, idperstart = 20, retry_max = 3), cores = c(e_mc = 1, m_mc = 1, nthreads = 1, blas = 1), verbose = TRUE)
Rprof(NULL)
summaryRprof()

library(profvis)
library(bench)
set.seed(123477)
profvis({
  group = sample(3, data[, length(unique(id))], replace = TRUE)
  fp = list(maxdf = 30, iter = 10, starts = 5, idperstart = 20, retry_max = 3)
  data[,group:=group[id]]
  f = trajectories(data, 3, group, fp = fp, cores = c(e_mc = 1, m_mc = 1, nthreads = 1, blas = 1), verbose = TRUE)
  if(!is.null( (er = xit_report(f, fp)) )) print(er)
})

