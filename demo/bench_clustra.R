library(clustra)

#' bench_clustra: Times several core use and data size configurations and
#' reports the results.
#' 
#' @param nid
#' Vector of data size configurations using number of id's
#' 
bench_clustra = function(nid = c(10000, 20000),
    fp = list(maxdf = 30, iter = 8, starts = 4, idperstart = 20, retry_max = 3),
                         c1 = 2, c2 = 2, c3 =  2, c4 = 2, 
                         verbose = FALSE) {
  rng_prev = RNGkind("L'Ecuyer-CMRG")
  for(n_id in nid) {
    data = gen_traj_data(n_id, m_obs = 25, s_range = c(-50, -10),
                         e_range = c(3*365, 10*365), reference = 100)
    results = data.frame(matrix(NA, nrow=length(nid)*c1*c2*c3*c4, ncol = 9))
    names(results) = c("n_id", "rep", "ec", "mc", "bc", "bl", 
                       "time", "xit", "all.equal")
    compare = NULL
    for(ec in seq(1, c1)) {
      for(mc in seq(1, c2)) {
        for(bc in seq(1, c3)) {
          for(bl in seq(1, c4)) {
            for(rep in 1:3) {
              openblasctl::openblas_set_num_threads(bl)
              set.seed(123473)
              cat("starting", n_id, rep, ec, mc, bc, bl, "...")
              time = system.time((
                cl = clustra(data, k = 3,
                             fp = list(maxdf = 30, iter = 10, starts = 5,
                                       idperstart = 20, retry_max = 3),
                             cores = c(e_mc = ec, m_mc = mc, nthreads = bc,
                                       blas = bl), verbose = verbose)
              ))
              er = paste(xit_report(cl), collapse = " ")

              if(is.null(compare)) {
                compare = cl$group
                agree = NA
              } else {
                agree = all.equal(compare, cl$group)
              }
              cat(time["elapsed"], agree, "\n")
            
              ir = (ec - 1)*c2*c3*c4 + (mc - 1)*c3*c4 + (bc - 1)*c4 + bl
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

bench_clustra(c(10000, 20000), 1, 1, 1, 4)