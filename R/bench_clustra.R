library(clustra)
bench_clustra = function(c1 = 4, c2 = 4, c3 =  4, c4 = 4, 
                         verbose = FALSE) {
  rng_prev = RNGkind("L'Ecuyer-CMRG")
  data = gen_traj_data(n_id = 1000, m_obs = 25, s_range = c(-50, -10),
                       e_range = c(3*365, 10*365), reference = 100)
  results = data.frame(matrix(NA, nrow=c1*c2*c3*c4, ncol = 5))
  names(results) = c("ec", "mc", "bc", "bl", "time", "all.equal")
  compare = NULL
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
                                   blas = bl), verbose = verbose)
          ))
          ## TODO add an all.equal on cl instances
          if(is.null(compare)) {
            compare = cl$group
          } else {
            agree = all.equal(compare, cl$group)
          }
          
          ir = (ec - 1)*c2*c3*c4 + (mc - 1)*c3*c4 + (bc - 1)*c4 + bl
          print((results[ir, ] = 
                   data.frame(ec = ec, mc = mc, bc = bc, bl = bl,
                              time = time["elapsed"], qll.equal = agree)))
        }
      }
    }
  }
  print(results)
}

bench_clustra(2, 2, 1, 1)