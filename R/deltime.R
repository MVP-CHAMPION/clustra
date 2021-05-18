#' Timing function
#'
#' @param ltime
#' Result of last call to deltime. 
#' @param text 
#' Text to display along with elapsed time.
#'
#' @return
#' "elapsed" component of current \code{\link{proc.time}}.
#' 
#' @export
deltime = function(ltime = proc.time()["elapsed"], text = NULL) {
  time = proc.time()["elapsed"]
  if(!is.null(text))
    cat(paste0(text, round(time - ltime, 1)))
  invisible(time)
}