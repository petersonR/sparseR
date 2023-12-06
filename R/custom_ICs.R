#' Custom IC functions for stepwise models
#'
#' @rdname custom_ics
#' @aliases EBIC RBIC RAIC
#' @param fit a fitted object
#' @param varnames names of variables
#' @param pen_info penalty information
#' @param gammafn What to use for gamma in formula
#' @param return_df should the deg. freedom be returned
#' @param ... additional args
#' @importFrom stats AIC BIC coef add.scope approx as.formula complete.cases
#'   deviance drop.scope extractAIC factor.scope formula glm lm.fit lm.wfit
#'   median model.frame model.matrix model.offset model.response model.weights
#'   nobs pf predict quantile terms update.formula weighted.residuals
#'   pchisq pf
#'
#' @return A vector of values for the criterion requested, and the degrees of
#'   freedom (appended to front of vector) if return_df == TRUE.
#' @export
#'
EBIC <- function(...){
  UseMethod("EBIC")
}

#' @rdname custom_ics
#' @export
EBIC.default <- function(fit, varnames, pen_info, gammafn = NULL, return_df = TRUE, ...) {

  if(!is.null(fit$model)) {
    n <- nrow(fit$model)
  } else if (!is.null(fit$nobs)) {
    n <- fit$nobs
  } else if (length(fit$residuals)) {
    n <- length(fit$residuals)
  }else stop("Must provide RBIC a sample size")


  P <- sum(pen_info$penalties[!duplicated(pen_info$vartypes)])

  beta <- fit$coefficients
  m <- length(beta) - 1

  if(length(gammafn)) {
    stopifnot(is.function(gammafn))
  } else {
    gammafn <- function(P, m, n) {
      (log(P / sqrt(n)) / log(P)) * (P >= sqrt(n))
    }
  }

  if(return_df) {
    return(c(sum(m), BIC(fit) + 2 * sum(gammafn(P, m, n) * lchoose(P, m))))
  }

  BIC(fit) + 2 * sum(gammafn(P, m, n) * lchoose(P, m))
}

#' @rdname custom_ics
#' @export
RBIC <- function(fit, ...){
  UseMethod("RBIC")
}

#' @rdname custom_ics
#' @export
RBIC.default <- function(fit, varnames, pen_info, gammafn = NULL, return_df = TRUE, ...) {

  if(!is.null(fit$model)) {
    n <- nrow(fit$model)
  } else if (!is.null(fit$nobs)) {
    n <- fit$nobs
  } else if (length(fit$residuals)) {
    n <- length(fit$residuals)
  }else stop("Must provide RBIC a sample size")


  P <- pen_info$penalties[!duplicated(pen_info$vartypes)]
  P_index <- lapply(
    pen_info$vartypes[!duplicated(pen_info$vartypes)],
    function(g) {
    varnames[pen_info$vartypes == g]
  })

  beta <- fit$coefficients
  m <- sapply(P_index, function(g) {
    sum(g %in% names(beta))
  })

  if(length(gammafn)) {
    stopifnot(is.function(gammafn))
  } else {
    gammafn <- function(P, m, n) {
      (log(P / sqrt(n)) / log(P)) * (P >= sqrt(n))
    }
  }

  if(return_df) {
    return(c(sum(m), BIC(fit) + 2 * sum(gammafn(P, m, n) * lchoose(P, m))))
  }

  BIC(fit) + 2 * sum(gammafn(P, m, n) * lchoose(P, m))
}

# Helpers (internal)
myBIC <- function(fit, ...) {
  m <- length(fit$coefficients) - 1
  c(m, BIC(fit))
}

myAIC <- function(fit, ...) {
  m <- length(fit$coefficients) - 1
  c(m, AIC(fit))
}

#' @rdname custom_ics
#' @export
RAIC <- function(fit, ...){
  UseMethod("RAIC")
}

#' @rdname custom_ics
#' @export
RAIC.default <- function(fit, varnames, pen_info, gammafn = NULL, return_df = TRUE, ...) {
  if(!is.null(fit$model)) {
    n <- nrow(fit$model)
  } else if (!is.null(fit$nobs)) {
    n <- fit$nobs
  } else if (length(fit$residuals)) {
    n <- length(fit$residuals)
  }else stop("Must provide RAIC a sample size")


  P <- pen_info$penalties[!duplicated(pen_info$vartypes)]
  P_index <- lapply(
    pen_info$vartypes[!duplicated(pen_info$vartypes)],
    function(g) {
      varnames[pen_info$vartypes == g]
    })

  beta <- fit$coefficients
  m <- sapply(P_index, function(g) {
    sum(g %in% names(beta))
  })

  Pstar <- pmin(P, n - 1)

  if(return_df) {
    return(c(sum(m), AIC(fit) + 2*m*((n-m)/(n-Pstar) - 1) - 2*((n-m)/(n-Pstar) - 1)^2))
  }

  AIC(fit) + 2*m*((n-m)/(n-Pstar) - 1) - 2*((n-m)/(n-Pstar) - 1)^2
}

#' @export
RBIC.sparseRBIC <- function(fit, ...) {
  RBIC(fit$fit, varnames = fit$pen_info$varnames, pen_info = fit$pen_info,
       return_df = FALSE,...)
}

#' @export
RAIC.sparseRBIC <- function(fit, ...) {
  RAIC(fit$fit, varnames = fit$pen_info$varnames, pen_info = fit$pen_info,
       return_df = FALSE, ...)
}

#' @export
EBIC.sparseRBIC <- function(fit, ...) {
  EBIC(fit$fit, varnames = fit$pen_info$varnames, pen_info = fit$pen_info,
       return_df = FALSE, ...)
}
