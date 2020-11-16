library(clustra)
bench_interactive = function(c1 = 4, c2 = 4, c3 =  4, c4 = 4) {
  set.seed(12765)
  rng_prev = RNGkind("L'Ecuyer-CMRG")
  data = gen_traj_data(n_id = 1000, m_obs = 25, s_range = c(-50, -10),
                       e_range = c(3*365, 10*365), reference = 100)
  results = data.frame(matrix(NA, nrow=c1*c2*c3*c4, ncol = 5))
  names(results) = c("ec", "mc", "bc", "bl", "time")
  for(ec in seq(1, c1)) {
    for(mc in seq(1, c2)) {
      for(bc in seq(1, c3)) {
        for(bl in seq(1, c4)) { # VECLIB_MAXIMUM_THREADS
          set.seed(123473)
          cat("starting", ec, mc, bc, bl, "...\n")
          time = system.time((
            cl = clustra(data, k = 3,
                         fp = list(maxdf = 30, iter = 10, starts = 5,
                                   idperstart = 20, retry_max = 3),
                         cores = c(e_mc = ec, m_mc = mc, nthreads = bc,
                                   blas = bl), verbose = TRUE)
          ))
          ## TODO add an all.equal on cl instances
          ir = (ec - 1)*c2*c3*c4 + (mc - 1)*c3*c4 + (bc - 1)*c4 + bl
          print((results[ir, ] = c(ec, mc, bc, bl, time["elapsed"])))
        }
      }
    }
  }
  results
}

## short test
sessionInfo()
set.seed(3)
dat1 <- mgcv::gamSim(1, n=25000, dist="normal", scale=20)
dat2 <- mgcv::gamSim(1, n=25000, dist="normal", scale=20)
bs <- "cr";k <- 12

b1 <- mgcv::bam(y ~ s(x0,bs=bs)+s(x1,bs=bs)+s(x2,bs=bs,k=k)+s(x3,bs=bs),data=dat1)
b2 <- mgcv::bam(y ~ s(x0,bs=bs)+s(x1,bs=bs)+s(x2,bs=bs,k=k)+s(x3,bs=bs),data=dat2)

dat <- list(dat1, dat2)
b <- lapply(dat, function(dat, bs, k) 
  mgcv::bam(y ~ s(x0,bs=bs)+s(x1,bs=bs)+s(x2,bs=bs,k=k)+s(x3,bs=bs),data=dat),
            bs=bs, k=k)
all.equal(b1$coefficients, b[[1]]$coefficients)
all.equal(b2$coefficients, b[[2]]$coefficients)

b <- parallel::mclapply(dat, function(dat, bs, k) 
  mgcv::bam(y ~ s(x0,bs=bs)+s(x1,bs=bs)+s(x2,bs=bs,k=k)+s(x3,bs=bs),data=dat),
  bs=bs, k=k, mc.cores = 2)

all.equal(b1$coefficients, b[[1]]$coefficients)
all.equal(b2$coefficients, b[[2]]$coefficients)
## works on Rhea!!
## could it be vecLib??
## try OpenBLAS on Rhea!
## try setting single core vecLib on the mac!



#' Function to benchmark various core combinations for components of clustra
#'
bench_batch = function() {

  ## start with single node optimization
  
  
  ## loop over things to vary    
  for(i in 0:max2) {
    ts = 0
    cores = 2^i
    openblasctl::openblas_set_num_threads(cores)
    for(rp in 1:reps) {
      row = as.vector(system.time(FUN(...)))[1:3]
      bench[[i*reps + rp]] = c(cores, rp, row)
      ts = ts + row[3]
    }
    cat("bench_cores: finished", cores, ts, "\n")
  }
  bench = data.frame(do.call(rbind, bench))
  names(bench) = c("cores", "rep", "user", "system", "elapsed")
  bench = gather(bench, key = "type", value = "seconds",
                 user, system, elapsed)  %>%
    mutate(cores = factor(cores))
  
  p = ggplot(bench, aes(cores, seconds, color = type, group = type)) +
    geom_point() + stat_summary(fun.y=mean, geom="line") +
    ggtitle(paste("Scaling", FUNname))
  
  pdf(paste0(FUNname, "_blas.pdf"), width = 7, height = 7)
  print(p)
  dev.off()
}

## TODO UNFINISHED!
stop("Demo is incomplete: don not run.")
## Edits a Slurm submission script to exercise a number of scaling options
## 

script = readLines("exec/trajectories.slurm")

## Set working directory
playdir = "~/clustra_play"
if(!dir.exists(playdir)) dir.create(playdir)

## Set or get clustra parameters
PL = clustra_par(playdir)

## Set cores for OpenBlas if available
if(cores$blas > 1) openblasctl::openblas_set_num_threads(cores$blas)

max_cores = list(e_mc = 4, m_mc = 4, bam_nthreads = 1, blas = 1)
cor_combo = do.call(expand.grid, lapply(max_cores, function(x) 1:x))
for(i in 1:nrow(cor_combo)) {
  
}
