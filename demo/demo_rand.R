## Example of running clustra on a single generated data set
## 

## create a directory to store a json control structure 
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
    gen_par = list(seed = 90, n_id = 10000, m_obs = 25,
                   e_range = c(365*3, 365*10), plots = FALSE),
    cores = list(e_mc = 4, m_mc = 4, bam_nthreads = 1, blas = 1),
    traj_par = list(seed = 79, maxdf = 30, iter = 8, 
                    starts = list(ns = 5, nid = 10),
                    k_vec = c(2, 3, 4), replicates =4)
  )
  write_json(PL, parname, pretty = TRUE)
}

## develop two examples: clustra and rand_clustra
## 
## the double loop below should be transferred to rand_clustra in trajectories
## !!!!!!!!!!!!!!!!!!!!!!!!!

set.seed(PL$gen_par$seed)
a0 = a = deltime()

data = gen_long_data(n_id = PL$gen_par$n_id, m_obs = PL$gen_par$m_obs,
                    e_range = PL$gen_par$e_range, plots = PL$gen_par$plots)
a = a_fit = deltime(a, paste0("\nData (", paste(dim(data), collapse = ","), ") generated"))

## TODO make this a function rand_clustra and make one iteration clustra. Export
##      both.  Then the clustra::: call can be removed! Finally, put a data gen
##      with a single clustra call in demo. Also a rand_clustra call in demo.
set.seed(PL$traj_par$seed)
cores = PL$cores
k_vec = PL$traj_par$k_vec
maxdf = PL$traj_par$maxdf
starts = PL$traj_par$starts
iter = PL$traj_par$iter
replicates = PL$traj_par$replicates
results = vector("list", replicates*length(k_vec))

for(j in 1:length(k_vec)) {
  k = k_vec[j]
  for(i in 1:replicates) {
    a_0 = deltime()
    
    f = clustra(data, k, starts, cores)
    results[[(j - 1)*replicates + i]] = list(k = as.integer(k), 
                                             rep = as.integer(i),
                                             deviance = f$deviance,
                                             group = f$group)
    data$group = as.factor(f$group[data$id])
    a = deltime()
    cat(k, i, "it =", f$iterations, "dev =", f$deviance, "err =", f$try_errors,
        "ch =", f$changes, "time =", a - a_0, "\n")
  }
}

## save object results and parameters
save(results, k_vec, maxdf, starts, iter, replicates, file = "results.Rdata")
a = deltime(a, "\nSaved results")

a = deltime(a_fit, "\nTotal Fit time")
a = deltime(a0, "\nTotal time")

## plot Rand Index evaluation
RandIndex_pairs = clustra:::allpair_RandIndex(results)
clustra:::rand_plot(RandIndex_pairs)
a = deltime(a, "\nRandIndex time")
