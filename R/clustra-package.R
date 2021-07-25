#' clustra-package
#' 
#' Clusters medical trajectories (unequally spaced and unequal lengths) aligned 
#' by an intervention time. Performs k-means clustering, where each mean is a 
#' thin plate spline fit to all points in a cluster. Distance is MSE across 
#' trajectory points to cluster spline. Provides silhouette plots and Adjusted
#' Rand Index evaluations of the number of clusters. Scales well to large data
#' with multicore parallelism available to speed computation.
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
