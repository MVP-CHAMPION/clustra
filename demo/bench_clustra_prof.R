## profile clustra
set.seed(123473)
data = gen_traj_data(80000, m_obs = 45, s_range = c(-50, -10),
                     e_range = c(3*365, 10*365), reference = 100)
Rprof()
cl = clustra(data, k = 3,
             fp = list(maxdf = 30, iter = 10, starts = 5,
                       idperstart = 20, retry_max = 3),
             cores = c(e_mc = 1, m_mc = 1, nthreads = 1,
                       blas = 1), verbose = FALSE)
Rprof(NULL)
summaryRprof()

