## Example of running a clustra ensemble with Rand index to choose k
## 
library(clustra)

## Set working directory
playdir = "~/clustra_play"
if(!dir.exists(playdir)) dir.create(playdir)

## Set or get clustra parameters
PL = clustra_par(playdir)

## Set seed for reproducibility
## TODO Check if reproducible when running parallel
set.seed(PL$gen_par$seed)
a0 = a = deltime()

data = gen_traj_data(PL)
a = a_fit = deltime(a, paste0("\nData (", paste(dim(data), collapse = ","), ") generated"))

## cluster the trajectories
rand_clustra(data, PL, verbose = TRUE)
