#' clustra-package
#' 
#' Clusters trajectories (unequally spaced and unequal length time series) on 
#' a common time axis. Clustering proceeds by an EM algorithm that iterates 
#' switching between fitting a thin plate spline (TPS) to combined responses 
#' within each cluster (M-step) and reassigning cluster membership based on the 
#' nearest fitted TPS (E-step). Initial cluster assignments are random or 
#' distant trajectories. The fitting is done with the *mgcv* package function 
#' *bam*, which scales well to very large data sets. Additional parallelism 
#' available via multicore on unix and mac platforms.
#' 
#' @name clustra-package
#' @docType package
#' @author George Ostrouchov, David Gagnon, Hanna Gerlovin
#' @keywords Package
#' 
#' # Import package operators
#' @importFrom data.table ":="
#' 
NULL
