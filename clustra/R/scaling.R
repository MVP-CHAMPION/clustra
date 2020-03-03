##
## Edits a Slurm submission script to exercise a number of scaling options
## 

script = readLines("trajectories.slurm")
parname = "trajectories.par" 
if(file.exists(parname)) {
  PL = read_json(parname, simplifyVector = TRUE)
} else {
  PL = list( # Default parameters
    gen_par = list(seed = 90, n_id = 20000, m_obs = 25,
                   e_range = c(365*3, 365*10), plots = FALSE),
    cores = list(e_mc = 4, m_mc = 4, bam_nthreads = 1, blas = 1),
    traj_par = list(seed = 79, maxdf = 50, iter = 20, starts = 8, idperstart = 20,
                    ng_vec = c(2, 3, 4, 5), replicates = 10)
  )
  write_json(PL, parname, pretty = TRUE)
}

max_cores = list(e_mc = 4, m_mc = 4, bam_nthreads = 1, blas = 1)
cor_combo = do.call(expand.grid, lapply(max_cores, function(x) 1:x))
for(i in 1:nrow(cor_combo)) {
  
}

