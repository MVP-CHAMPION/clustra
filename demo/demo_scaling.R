## TODO UNFINISHED!
## Edits a Slurm submission script to exercise a number of scaling options
## 

script = readLines("trajectories.slurm")

PL = clustra_par(parname = "clustra.par", playdir = "~/clustra_play")

max_cores = list(e_mc = 4, m_mc = 4, bam_nthreads = 1, blas = 1)
cor_combo = do.call(expand.grid, lapply(max_cores, function(x) 1:x))
for(i in 1:nrow(cor_combo)) {
  
}

