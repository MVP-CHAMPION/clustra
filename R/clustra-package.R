#' clustra-package
#' 
#' Performs k-means clustering, where each mean is a thin plate spline
#' fit to all points of trajectories in cluster. 
#' 
#' @import ggplot2
#' @importFrom graphics abline axis box image layout par plot
#' @importFrom jsonlite read_json write_json
#' @importFrom openblasctl openblas_set_num_threads
#' @importFrom parallel mclapply
#' @importFrom grDevices dev.off pdf colorRampPalette
#' @importFrom stats dist predict rnorm rpois runif time
#' @importFrom utils write.table
#' 
#' @name clustra-package
#' @docType package
#' @author George Ostrouchov, David Gagnon, Hanna Gerlovin
#' @keywords Package
NULL
