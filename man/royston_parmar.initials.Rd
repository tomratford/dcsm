% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/royston_parmar.R
\name{royston_parmar.initials}
\alias{royston_parmar.initials}
\title{Get a set of (hopefully) feasible inital values for a Natural Cubic Spline based model}
\usage{
royston_parmar.initials(data, k01, k02, k12)
}
\arguments{
\item{data}{Data set to model from, with columns `L`, `R`, `V`, `Delta1`, `Delta2`}

\item{k01}{Number of interior knots in the 0 -> 1 transition}

\item{k02}{Number of interior knots in the 0 -> 2 transition}

\item{k12}{Number of interior knots in the 1 -> 2 transition}
}
\value{
A list of parameters, to be fed to \link{royston_parmar.ll}
}
\description{
Get a set of (hopefully) feasible inital values for a Natural Cubic Spline based model
}
