library(clustra)

#' bench_clustra: Times several core use and data size configurations and
#' reports the results.
#' 
#' OpenBLAS lib gives ~30% improvement on 80K-45 but no multicore improvement.
#' 
#' @param nid
#' Vector of data size configurations using number of id's
#' 
bench_clustra_core = function(nid,
                         e_mc = 1:2, m_mc = 1:2, nthreads =  1:2, blas = 1, reps = 1:2,
    fp = list(maxdf = 30, iter = 8, starts = 4, idperstart = 20, retry_max = 3),
                         verbose = FALSE) {
  rng_prev = RNGkind("L'Ecuyer-CMRG")
  for(n_id in nid) {
    set.seed(123473)
    data = gen_traj_data(n_id, m_obs = 45, s_range = c(-50, -10),
                         e_range = c(3*365, 10*365), reference = 100)
    rows = length(nid)*
      length(e_mc)*length(m_mc)*length(nthreads)*length(blas)*
      length(reps)
    results = data.frame(matrix(NA, nrow = rows, ncol = 9))
    names(results) = c("n_id", "rep", "ec", "mc", "bc", "bl", 
                       "time", "xit", "all.equal")
    compare = NULL
    for(ec in e_mc) {
      for(mc in m_mc) {
        for(bc in nthreads) {
          for(bl in blas) {
            for(rep in reps) {
              openblasctl::openblas_set_num_threads(bl)
              set.seed(1234737)
              cat("starting", n_id, rep, ec, mc, bc, bl, "...")
              time = system.time((
                cl = clustra(data, k = 3,
                             fp = list(maxdf = 30, iter = 10, starts = 5,
                                       idperstart = 20, retry_max = 3),
                             cores = c(e_mc = ec, m_mc = mc, nthreads = bc,
                                       blas = bl), verbose = verbose)
              ))
              er = paste(xit_report(cl, fp), collapse = " ")

              if(is.null(compare)) {
                compare = cl$group
                agree = NA
              } else {
                agree = all.equal(compare, cl$group)
              }
              cat(er, time["elapsed"], agree, "\n")
            
              ir = (((((n_id - 1)*length(e_mc) + (ec - 1))*length(m_mc) + 
                (mc - 1))*length(nthreads) + (bc - 1))*length(blas) +
                (bl - 1))*length(reps) + rep
              results[ir, ] = data.frame(n_id = n_id, rep = rep,
                                         ec = ec, mc = mc, bc = bc, bl = bl,
                                         time = time["elapsed"], xit = er,
                                         qll.equal = agree)
            }
          }
        }
      }
    }
  }
  print(results)
}

# bench_clustra(c(80000, 160000), 1:6, 1:6, 1, 1, 1:2)
# bench_clustra(c(80000, 160000, 320000), 2:6, 2:6, 1, 1, 1:2)
