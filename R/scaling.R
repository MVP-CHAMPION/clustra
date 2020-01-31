cores = list(e_mc = 4, m_mc = 4, bam_nthreads = 1, blas = 1)
cores_fun = function(max_cores, cores = NULL) {
  ## will run all combinations of up to max_cores
  ## Use:
  ## cores = NULL 
  ## while( (cores = cores_fun(max_cores, cores)) ) {
  ##   function body
  ## }
  ##
  if(is.null(cores)) {
    cores = list(e_mc = 1, m_mc = 1, bam_nthreads = 1, blas = 1)
  } else {
    if(cores$e_mc < max_cores$e_mc) {
      cores$e_mc = cores$e_mc + 1
    } else {
      cores$e_mc = 1
      if(cores$m_mc < max_cores$m_mc) {
        cores$m_mc = cores$m_mc + 1
      } else {
        cores$m_mc = 1
        if(cores$bam_nthreads < max_cores$bam_nthreads) {
          cores$bam_nthreads = cores$bam_nthreads + 1
        } else {
          cores$bam_nthreads = 1
          if(cores$blas < max_cores$blas) {
            cores$blas = cores$blas + 1
          } else {
            ## done!
            cores = FALSE
          }
        }
      }
    }
  }
  cores
}
