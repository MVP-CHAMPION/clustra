
source("R/generate.R")
source("R/scaling.R")
source("R/evaluate.R")
source("R/trajectories.R")
set.seed(PL$gen_par$seed)

dat = gen_long_data(n_id = PL$gen_par$n_id, m_obs = PL$gen_par$m_obs,
                    e_range = PL$gen_par$e_range, plots = PL$gen_par$plots)
a = a_fit = deltime(a, paste0("\nData (", paste(dim(dat), collapse = ","), ") generated"))

set.seed(PL$traj_par$seed)
cores = PL$cores
ng_vec = PL$traj_par$ng_vec
maxdf = PL$traj_par$maxdf
starts = PL$traj_par$starts
iter = PL$traj_par$iter
start_nid = PL$traj_par$idperstart*length(ng_vec)
replicates = PL$traj_par$replicates
results = vector("list", replicates*length(ng_vec))
for(j in 1:length(ng_vec)) {
  ng = ng_vec[j]
  for(i in 1:replicates) {
    a = a_i = deltime(a)
    
    ## get initial group assignment
    group = start_groups(dat, ng, starts = starts, start_nid = start_nid,
                         cores = cores, maxdf = maxdf)
    
    dat$id = as.numeric(factor(dat$id)) # since id are not sequential
    dat$group = group[dat$id] # expand to all responses
    a = deltime(a, " Starts time")
    
    ## now dat includes initial group assignment
    f = trajectories(dat = dat, ng, group, iter = iter, maxdf = maxdf, plot = FALSE,
                     cores = cores)
    results[[(j - 1)*replicates + i]] = list(ng = as.integer(ng), 
                                             rep = as.integer(i),
                                             deviance = f$deviance,
                                             group = f$group)
    dat$group = as.factor(f$group[dat$id])
    a = deltime(a_i, " Replicate time")
  }
}

## save object results and parameters
save(results, ng_vec, maxdf, starts, iter, start_nid, replicates, file = "results.Rdata")
a = deltime(a, "\nSaved results")

a = deltime(a_fit, "\nTotal Fit time")
a = deltime(a0, "\nTotal time")

## plot Rand Index evaluation
RandIndex_pairs = allpair_RandIndex(results)
rand_plot(RandIndex_pairs)
a = deltime(a, "\nRandIndex time")
