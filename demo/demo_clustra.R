## Example of running clustra on a single generated data set
## 
library(clustra)
PL = clustra_par()

## Set seed for reproducibility
## TODO Check if reproducible when running parallel
set.seed(PL$gen_par$seed)
a0 = a = deltime()

## Generate a data set
data = gen_traj_data(n_id = PL$gen_par$n_id, lambda_obs = PL$gen_par$lambda_obs,
                    first = PL$gen_par$first, last = PL$gen_par$last,
                    plots = PL$gen_par$plots, wfile = "clustra_gen_data.csv")
a = a_fit = deltime(a, paste0("\nData (", paste(dim(data), collapse = ","), ") generated"))

## cluster the trajectories
cl = clustra(data, 3, PL, verbose = TRUE, plot = TRUE)

