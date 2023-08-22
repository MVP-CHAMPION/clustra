#' Simulated blood pressure data
#' 
#' Contains 80,000 individuals, each with an average of about 17 observations
#' in 5 clusters with scatter. Generated from a 5-cluster thin spline model of 
#' actual blood pressures collected from roughly the same number of individuals
#' at U.S. Department of Veterans Affairs facilities in connection with the
#' CHAMPION project. Each cluster-mean generated individual has a random number
#' of observations at random times with one observation at intervention time 0,
#' and with added standard normal error.
#' 
#' @format ## `bp`
#' A "data.table" and "data.frame" with 1,353,910 rows and 4 columns:
#' \describe{
#'   \item{id}{An integer in 1:80000.}
#'   \item{true_group}{An integer in 1:5.}
#'   \item{time}{An integer between -365 and 730, giving observation day with
#'   reference to an intervention at time 0.}
#'   \item{response}{The systolic blood pressure on that day.}
#' }
#' 
"bp"