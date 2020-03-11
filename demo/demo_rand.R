## Example of running a clustra ensemble with Rand index to choose k
## 

## Create a directory to store a json control structure 
playdir = "~/clustra_play"
if(!dir.exists(playdir)) dir.create(playdir)
setwd(playdir)

library(clustra)
library(jsonlite)
parname = "trajectories.par" 
if(file.exists(parname)) {
  PL = read_json(parname, simplifyVector = TRUE)
} else {
  PL = list( # Default parameters
    gen_par = list(seed = 90, n_id = 10000, lambda_obs = 25, first = c(-50, -10),
                   last = c(365*3, 365*10), plots = FALSE),
    cores = list(e_mc = 4, m_mc = 4, bam_nthreads = 1, blas = 1),
    traj_par = list(seed = 79, maxdf = 30, iter = 8, 
                    starts = list(ns = 5, nid = 10),
                    k_vec = c(2, 3, 4), replicates =4)
  )
  write_json(PL, parname, pretty = TRUE)
}

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
