#' clustra-package
#' 
#' Performs k-means clustering, where each mean is a thin plate spline fit to
#' all points from trajectories in cluster. Distance is the mean squared error
#' or maximum error of spline to trajectory points.
#' 
#' 
#' @name clustra-package
#' @docType package
#' @author George Ostrouchov, David Gagnon, Hanna Gerlovin
#' @keywords Package
#' 
#' # Import package operators
#' @importFrom data.table ":="

# Make sure data.table knows we know we're using it
.datatable.aware = TRUE

#' 
NULL
