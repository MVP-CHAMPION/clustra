% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trajectories.R
\name{pred_g}
\alias{pred_g}
\title{Function to predict for new data based on fitted gam object.}
\usage{
pred_g(tps, newdata)
}
\arguments{
\item{tps}{Output structure of \code{\link[mgcv]{bam}}.}

\item{newdata}{See \code{\link{clustra}} description of data parameter.}
}
\value{
A numeric vector of predicted values corresponding to rows of newdata.
If gam object is NULL, NULL is returned instead.
}
\description{
Function to predict for new data based on fitted gam object.
}
