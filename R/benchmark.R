#' Function to benchmark various core combinations for components of clustra
#'
#' @param FUN Function to benchmark.
#' @param max2 Largest power of 2 in cores (e.g. max2 = 3 implies 1, 2, 4, and 8
#'  cores will be tested. Core control is via openblasctl package)
#' @param reps Number of replicates for each core setting
#' @export
bench_cores = function(FUN, ..., max2 = 5, reps = 4) {
  require(openblasctl)
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  FUNname = substitute(FUN)
  FUN = match.fun(FUN)
  
  bench = vector("list", max2*reps)
  
  for(i in 0:max2) {
    ts = 0
    cores = 2^i
    openblas_set_num_threads(cores)
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

