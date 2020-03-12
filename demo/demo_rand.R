## Example of running a clustra ensemble with Rand index to choose k
## 

library(clustra)
library(jsonlite)

PL = clustra_par()

## Set seed for reproducibility
## TODO Check if reproducible when running parallel
set.seed(PL$gen_par$seed)
a0 = a = deltime()

data = gen_traj_data(n_id = PL$gen_par$n_id, lambda_obs = PL$gen_par$lambda_obs,
                    first = PL$gen_par$first, last = PL$gen_par$last,
                    plots = PL$gen_par$plots)
a = a_fit = deltime(a, paste0("\nData (", paste(dim(data), collapse = ","), ") generated"))

## cluster the trajectories
rand_clustra(data, PL, verbose = TRUE)
