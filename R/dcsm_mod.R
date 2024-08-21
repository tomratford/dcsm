#' Print the treatment effect estimates for a dual-censored model
#'
#' @param x A \code{dcsm_mod} object
#' @param ... Not used
#'
#' @return Nothing, prints three estimates of the treatment effect for each transition
#' @export
print.dcsm_mod <- function(x, ...) {
  thetas <- x$par[c("theta01","theta02","theta12")]
  cat(paste(names(thetas),collapse = "\t"))
  cat("\n")
  cat(paste(format(thetas,digits=3),collapse="\t"))
}

#' Plot the survivor curve for a dual-censored model
#'
#' @param x A \code{dcsm_mod} object
#' @param treat_col The colour of the treatment effect line, set to \code{NULL} to only plot the baseline survivor function.
#' @param add Add this to a existing plot?
#' @param ... Passed to \link{curve}
#'
#' @return Nothing, plots to output
#' @export
plot.dcsm_mod <- function(x,treat_col="blue",add=F,...) {
  if (attr(x,"dist") == "weibull") {
    pars <- as.list(x$par)
  } else {
    if (exists(attr(x,"initials"))) {
      pars <- make_pars2(x$par, attr(x,"initials"))
    } else {
      stop("No initials found, running this will brick your R session")
    }
  }
  fns <- switch(attr(x,"dist"),
                "weibull" = do.call(weib.fnBuilder,pars),
                "joly" = do.call(joly.fnBuilder,pars),
                "royston_parmar" = do.call(royston_parmar.fnBuilder,pars))
  curve(fns$P00(rep(0,length(x)),x,rep(0,length(x))),add=add,...)
  if (!is.null(treat_col))
    curve(fns$P00(rep(0,length(x)),x,rep(1,length(x))),add=T,col=treat_col,...)
}
