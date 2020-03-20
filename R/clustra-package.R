#' clustra-package
#' 
#' Performs k-means clustering, where each mean is a thin plate spline
#' fit to all points of trajectories in cluster. 
#' 
#' @importFrom dplyr mutate
#' @import ggplot2
#' @importFrom graphics abline axis box image layout par plot
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom jsonlite read_json write_json
#' @importFrom magrittr %>%
#' @importFrom openblasctl openblas_set_num_threads
#' @importFrom parallel mclapply
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats dist predict rnorm rpois runif time
#' @importFrom tidyr gather
#' @importFrom utils write.table
#' 
#' @name clustra-package
#' @docType package
#' @author George Ostrouchov, David Gagnon, Hanna Gerlovin
#' @keywords Package
NULL
