
#' Components of `gen` list in .clustra_env environment
#' 
#' @param seed Seed to use in data generation.
#' @param n_id Number of id to generate.
#' @param m_obs Mean number of observation per `id` using
#' @param s_range A vector of length 2, giving the min and max limits of
#'   uniformly generated start observation time.
#' @param e_range A vector of length 2, giving the min and max limits of
#'   uniformly generated end observation time.
#' @param plots Logical to indicate if a sample of `id`s trajectories are
#'   plotted.
#' @rdname env_gen
clustra_gen = function(
  seed = 90,
  n_id = 20000,
  m_obs = 25,
  s_range = c(-50, -10),
  e_range = c(1095, 3650),
  plots = FALSE
) {
  list(
    seed = seed,
    n_id = n_id,
    m_obs = m_obs,
    s_range = s_range,
    e_range = e_range,
    plots = plots
  )
}

#' Components of `cor` list in .clustra_env environment
#' 
#' @param e_mc Cores to use by mclapply of e_mc (expectation over k)
#' @param m_mc Cores to use by mclapply of m_mc (maximization across k)
#' @param bam_nthreads Threads to use by mgcv::bam()
#' @param blas Cores to use by OpenBLAS
#' @rdname env_cor
clustra_cor = function(
  e_mc = 1, 
  m_mc = 1,
  bam_nthreads = 1,
  blas = 1
) {
  list(
    e_mc = e_mc,
    m_mc = m_mc,
    bam_nthreads = bam_nthreads,
    blas = blas
  )
}

#' Components of `clu` list in .clustra_env environment
#' 
#' @param seed Cores to use by mclapply of e_mc (expectation over k)
#' @param maxdf Basis dimension of smooth term (mgcv::s(), parameter k).
#' @param iter Maximum number of iterations.
#' @param starts Number of random starts.
#' @param idperstart Number of `id`s to sample for random starts.
#' @rdname env_cor
#' @seealso [mgcv::s()]
clustra_clu = function(
  seed = 79, # starting seed for cluster starts
  maxdf = 50, # 
  iter = 8, # mgcv::bam maximum iterations
  starts = 4,
  idperstart = 20
) {
  list(
    seed = seed,
    maxdf = maxdf,
    iter =iter,
    starts = starts,
    idperstart = idperstart
  )
}

clustra_ran = function(
  ng_vec = c(2, 3, 4, 5),
  replicates = 10
)
{
  list(
    ng_vec = ng_vec,
    replicates = replicates
  )  
}

#' Function to get or set clustra environment variables.
#' 
#' @param ... Any number of string constants that ask for or assign various
#' parameters in .clustra_env environment. For example, "cor$e_mc" will
#' return the cor$e_mc parameter and "cor$e_mc = 8" will assign 8 to that
#' parameter.
#' @details 
#' Several parameters can be set but only one parameter returned. All
#'  parameters after the first returned parameter are ignored. An empty
#'  parameter list will print all parameters in the environment.
#'  Note that limited checking for syntax is performed so be careful!
#' @export
clustra_env = function(...) {
  fun = "clustra_env"
  env = .GlobalEnv$.clustra_env
  enc = ".clustra_env"
  args = list(...)
  for(a in args) {
    if(!is.character(a)) {
      cat(fun, ": expecting character string, but got:\n")
      print(a)
      return(invisible())
    }
    ax = gsub(" ", "", a, fixed = TRUE)
    alr = unlist(strsplit(ax, split = "=", fixed = TRUE)) # separate assignment
    if(length(alr) > 2) {
      cat(fun, ": too many = in", a, "\n")
      return(invisible())
    }
    al = unlist(strsplit(alr[1], split = "$", fixed = TRUE)) # separate reference
    if(length(al) > 2) {
      cat(fun, ": too many $ in", a, "\n")
      return(invisible())
    }
    if(!exists(al[1], envir = env)) { # object in env?
      cat(fun, ": did not find object", al[1], "in", enc, "\n")
      return(invisible())
    }
    if(length(al) > 1 & is.null(env[[al[1]]][[al[2]]])) { # component in object?
      cat(fun, ": component", al[2], "not in", al[1], "of", enc, "\n")
      return(invisible())
    }
    
    if(length(alr) == 1) {
      ## return requested parameter and ignore rest of parameters, if any
      return(eval(parse(text = a), envir = env))
      } else if(length(alr) == 2) {
        ## silently modify environment parameter
        eval(parse(text = a), env)
      } else {
        cat(fun, ": something went wrong with parameter", a, "\n")
        return(invisible())
      }
  }
  
  ## an empty parameter list prints all parameters in environment
  if(length(a) == 0)
    print(sapply(ls(env), function(x) get(x, envir = env)))
    
  invisible()
}

### Initialize clustra options
.clustra_init = function(envir = .GlobalEnv){
  if(!exists(".clustra_env", envir = envir)){
    envir$.clustra_env = new.env()
  } 
  envir$.clustra_env$gen = clustra_gen()
  envir$.clustra_env$cor = clustra_cor()
  envir$.clustra_env$clu = clustra_clu()
  envir$.clustra_env$ran = clustra_ran()
  
  invisible()
}

.onLoad = function(libname, pkgname) {
  .clustra_init()
}
.onUnload = function(libpath) {
  rm(.clustra_env)
}
