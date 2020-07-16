## Example of running clustra on a single generated data set
##
library(clustra)

## Set working directory
playdir = "~/clustra_play"
if(!dir.exists(playdir)) dir.create(playdir)

## Generate data
set.seed(clustra_env("gen$seed"))
a0 = a = deltime()
data = gen_traj_data()

## cluster the trajectories
clustra_env("clu$iter = 10", "clu$seed = 123473",
            "cor$e_mc = 1", "cor$m_mc = 1")
set.seed(clustra_env("clu$seed"))
cl = clustra(data, 3, verbose = TRUE)
