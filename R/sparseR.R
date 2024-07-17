#' Fit a ranked-sparsity model with regularized regression
#'
#' @param formula Names of the terms
#' @param data Data
#' @param family The family of the model
#' @param penalty What penalty should be used (lasso, MCP, or SCAD)
#' @param alpha The mix of L1 penalty (lower values introduce more L2 ridge
#'   penalty)
#' @param ncvgamma The tuning parameter for ncvreg (for MCP or SCAD)
#' @param lambda.min The minimum value to be used for lambda (as ratio of max,
#'   see ?ncvreg)
#' @param k The maximum order of interactions to consider (default: 1; all
#'   pairwise)
#' @param poly The maximum order of polynomials to consider (default: 2)
#' @param gamma The degree of extremity of sparsity rankings (see details)
#' @param cumulative_k Should penalties be increased cumulatively as order
#'   interaction increases?
#' @param cumulative_poly Should penalties be increased cumulatively as order
#'   polynomial increases?
#' @param pool Should interactions of order k and polynomials of order k+1 be
#'   pooled together for calculating the penalty?
#' @param ia_formula formula to be passed to step_interact (for interactions,
#'   see details)
#' @param pre_process Should the data be preprocessed (if FALSE, must provide
#'   model_matrix)
#' @param model_matrix A data frame or matrix specifying the full model matrix
#'   (used if !pre_process)
#' @param y A vector of responses (used if !pre_process)
#' @param poly_prefix If model_matrix is specified, what is the prefix for
#'   polynomial terms?
#' @param int_sep If model_matrix is specified, what is the separator for
#'   interaction terms?
#' @param pre_proc_opts List of preprocessing steps (see details)
#' @param filter The type of filter applied to main effects + interactions
#' @param extra_opts A list of options for all preprocess steps (see details)
#' @param ... Additional arguments (passed to fitting function)
#'
#' @details
#'
#' Selecting `gamma`: higher values of gamma will penalize "group" size more. By
#' default, this is set to 0.5, which yields equal contribution of prior
#' information across orders of interactions/polynomials (this is a good
#' default for most settings).
#'
#' Additionally, setting `cumulative_poly` or `cumulative_k` to `TRUE` increases
#' the penalty cumulatively based on the order of either polynomial or
#' interaction.
#'
#' The options that can be passed to `pre_proc_opts` are: - knnImpute (should
#' missing data be imputed?) - scale (should data be standardized)? - center
#' (should data be centered to the mean or another value?) - otherbin (should
#' factors with low prevalence be combined?) - none (should no preprocessing be
#' done? can also specify a null object)
#'
#' The options that can be passed to `extra_opts` are:
#' - centers (named numeric vector which denotes where each covariate should be
#'   centered)
#' - center_fn (alternatively, a function can be specified to calculate center
#'   such as `min`
#' or `median`)
#' - freq_cut, unique_cut (see ?step_nzv; these get used by the filtering steps)
#' - neighbors (the number of neighbors for knnImpute)
#' - one_hot (see ?step_dummy), this defaults to cell-means coding which can be
#' done in regularized regression (change at your own risk)
#' - raw (should polynomials not be orthogonal? defaults to true because variables are
#' centered and scaled already by this point by default)
#'
#' \code{ia_formula} will by default interact all variables with each other up
#' to order k. If specified, ia_formula will be passed as the `terms` argument
#' to \code{recipes::step_interact}, so the help documentation for that function
#' can be investigated for further assistance in specifying specific
#' interactions.
#'
#' @importFrom ncvreg cv.ncvsurv cv.ncvreg
#' @importFrom rlang .data
#' @importFrom dplyr everything group_by summarize n ungroup mutate
#'
#' @md
#'
#' @return an object of class `sparseR` containing the following:
#'
#' \item{fit}{the fit object returned by `ncvreg`}
#' \item{srprep}{a `recipes` object used to prep the data}
#' \item{pen_factors}{the factor multiple on penalties for ranked sparsity}
#' \item{results}{all coefficients and penalty factors at minimum CV lambda}
#' \item{results_summary}{a tibble of summary results at minimum CV lambda}
#' \item{results1se}{all coefficients and penalty factors at lambda_1se}
#' \item{results1se_summary}{a tibble of summary results at lambda_1se}
#' \item{data}{the (unprocessed) data}
#' \item{family}{the family argument (for non-normal, eg. poisson)}
#' \item{info}{a list containing meta-info about the procedure}
#'
#' @references For fitting functionality, the `ncvreg` package is used; see
#' Breheny, P. and Huang, J. (2011) Coordinate descent algorithms for nonconvex
#' penalized regression, with applications to biological feature selection. Ann.
#' Appl. Statist., 5: 232-253.
#'
#' @export

sparseR <- function(formula, data, family = c("gaussian", "binomial", "poisson", "coxph"),
                    penalty = c("lasso", "MCP", "SCAD"), alpha = 1, ncvgamma = 3,
                    lambda.min = .005,
                    k = 1, poly = 2, gamma = .5, cumulative_k = FALSE,
                    cumulative_poly = TRUE, pool = FALSE,
                    ia_formula = NULL,
                    pre_process = TRUE, model_matrix = NULL, y = NULL,
                    poly_prefix = "_poly_", int_sep = "\\:",
                    pre_proc_opts = c("knnImpute", "scale", "center", "otherbin", "none"),
                    filter = c("nzv", "zv"), extra_opts = list(), ...) {

  ## Check and validate arguments
  if(is.null(pre_proc_opts))
    pre_proc_opts <- "none"
  pre_proc_opts <- match.arg(pre_proc_opts, several.ok = TRUE)
  filter <- match.arg(filter)
  penalty <- match.arg(penalty)
  family <- match.arg(family)

  if(family == "coxph") {

    fit_fn <- ncvreg::cv.ncvsurv
  } else {
    fit_fn <- ncvreg::cv.ncvreg
    if(!is.null(y) && ("Surv" %in% class(y)))
      stop("Detected survival outcome, do you mean to set family = 'coxph'?")
  }


  if(pre_process & !is.null(model_matrix))
    stop("If specifying model matrix, pre_process must be set to FALSE")

  if(pre_process) {

    ## Create and prep model matrix function (preprocessing)
    srprep <- sparseR_prep(formula, data, k, poly, pre_proc_opts,
                           filter, extra_opts, ia_formula = ia_formula,
                           family = family)

    ## "Bake" recipe using preprocess steps (this can be done on test data too later)
    X <- bake(srprep, data, everything(), - all_outcomes())

    # Extract outcome
    if(!is.null(srprep$outcome)) {
      y <- srprep$outcome
    } else {
      stop("formula must have outcome for sparseR")
    }

  }  else {
    srprep <- NULL
    if(is.null(model_matrix))
      stop("If pre_process is FALSE, model_matrix must be specified")

    if(is.null(y))
      stop("If pre_process is FALSE, the outcome y must be specified")

    stopifnot(is.data.frame(as.data.frame(model_matrix)))
    X <- as.data.frame(model_matrix)
    cc <- complete.cases(X) & complete.cases(y)

    if(any(!cc))
      warning("Omitting ", sum(!cc), " missing observations")

    X <- X[cc,]
    y <- y[cc]

    data <- data.frame(y, X)
    names(data)[1] <- deparse(substitute(y))[[1]]
  }

  ## Get the correct penalty scaling for each covariate
  pen_info <- get_penalties(
    names(X), poly, pool = pool, cumulative_poly = cumulative_poly,
    cumulative_k = cumulative_k, gamma = gamma, poly_prefix = poly_prefix,
    int_sep = int_sep
  )

  pn <- pen_info$penalties

  ## Run regularized regression model
  fit <- fit_fn(as.matrix(X), y, family = family, penalty = penalty,
                   penalty.factor = pn, alpha = alpha, lambda.min = lambda.min,
                   gamma = ncvgamma, ...)

  ## Get results @ CV MIN
  coefs <- predict(fit, type = "coef")[,1]
  if(family != "coxph")
    coefs <- coefs[-1]

  results <- data.frame(
    "Vartype" = c(pen_info$vartypes),
    coef = coefs,
    penalty = c(pn)
  )

  results_summary <- results %>%
    group_by(.data$Vartype) %>%
    summarize(Total = n(), Selected = sum(.data$coef != 0),
              Saturation = .data$Selected/.data$Total,
              Penalty = mean(.data$penalty, na.rm = TRUE))%>%
    ungroup() %>%
    mutate(Vartype = as.character(.data$Vartype))

  ## Get results @ CV 1se
  cv1se <- fit$cve + fit$cvse
  coefs <- predict(fit, type = "coef", which = which(fit$cve < min(cv1se))[[1]])[,1]
  if(family != "coxph")
    coefs <- coefs[-1]

  results1se <- data.frame(
    "Vartype" = c(pen_info$vartypes),
    coef = coefs,
    penalty = c(pn)
  )

  results1se_summary <- results1se %>%
    group_by(.data$Vartype) %>%
    summarize(Total = n(), Selected = sum(.data$coef != 0),
              Saturation = .data$Selected/.data$Total,
              Penalty = mean(.data$penalty, na.rm = TRUE))%>%
    ungroup() %>%
    mutate(Vartype = as.character(.data$Vartype))

  info <- list(k = k, poly = poly, cumulative_k = cumulative_k,
               cumulative_poly = cumulative_poly, pool = pool,
               pre_process = pre_process, model_matrix = model_matrix, y = y,
               poly_prefix = poly_prefix, int_sep = int_sep)


  ## Return relevant results
  res <- list(fit = fit, srprep = srprep, pen_factors = pn,
              results = results, results_summary = results_summary,
              results1se = results1se, results1se_summary = results1se_summary,
              data = data, family = family, info = info)
  class(res) <- "sparseR"
  res
}

#' Print sparseR object
#'
#' @param x a sparseR object
#' @param prep Should the SR set-up information be printed as well?
#' @param ... additional arguments passed to print.ncvreg
#'
#' @method print sparseR
#'
#' @return returns x invisibly
#'
#' @export
print.sparseR <- function(x, prep = FALSE, ...) {
  if(prep){
    cat("Model matrix setup information:\n\n")
    print(x$srprep)
  }
  cat("\nModel summary @ min CV:\n")
  cat("-----------------------------------------------------\n")
  alt_print(suppressMessages(summary(x$fit)))
  cat("\n  SR information:\n")

  print(as.data.frame(x$results_summary), row.names = FALSE, digits = 3)

  cat("\n\nModel summary @ CV1se:\n")
  cat("-----------------------------------------------------\n")
  alt_print(summary(x$fit), which = which(x$fit$cve < min(x$fit$cve + x$fit$cvse))[1])
  cat("\n  SR information:\n")
  print(as.data.frame(x$results1se_summary), row.names = FALSE, digits = 3)
  return(invisible(x))
}

## Alternate printing method for summary.cv.ncvreg that allows for different lambda val
alt_print <- function (x, digits, which, ...) {

  if(missing(which))
    which <- x$min
  digits <- if (missing(digits))
    digits <- c(2, 4, 2, 2, 3)
  else rep(digits, length.out = 5)
  cat("  ", x$penalty, "-penalized ", x$model, " regression with n=",
      x$n, ", p=", x$p, "\n", sep = "")
  cat("  (At lambda=", formatC(x$lambda[which], digits[2], format = "f"), "):\n", sep = "")
  cat("    Nonzero coefficients: ", x$nvars[which], "\n", sep = "")
  cat("    Cross-validation error (deviance): ",
      formatC(x$cve[which], digits[1], format = "f"), "\n", sep = "")
  cat("    R-squared: ", formatC(x$r.squared[which], digits[3], format = "f"),
      "\n", sep = "")
  cat("    Signal-to-noise ratio: ",
      formatC(x$snr[which], digits[4], format = "f"), "\n", sep = "")
  if (x$model == "logistic")
    cat("    Prediction error: ", formatC(x$pe[which], digits[5], format = "f"),
        "\n", sep = "")
  if (x$model == "linear")
    cat("    Scale estimate (sigma): ",
        formatC(sqrt(x$cve[which]),  digits[5], format = "f"), "\n", sep = "")
}


#' Predict coefficients or responses for sparseR object
#'
#' @param object sparseR object
#' @param newdata new data on which to make predictions
#' @param lambda a particular value of lambda to predict with
#' @param at a "smart" guess to use for lambda
#' @param ... additional arguments passed to predict.ncvreg
#'
#' @method predict sparseR
#'
#' @return predicted outcomes for `newdata` (or coefficients)
#'   at specified (or smart) lambda value
#'
#' @export
predict.sparseR <- function(object, newdata, lambda, at = c("cvmin", "cv1se"), ...) {

  at <- match.arg(at)

  args <- list(...)

  ## If newdata is not specified...
  if(missing(newdata)) {
    if(!is.null(object$srprep)) {
      X <- as.matrix(bake(object$srprep, object$data, everything(), -all_outcomes()))
    } else {
      X <- NULL
    }
  } else {
    if(!is.null(object$srprep)) {
      outcome_name <- object$srprep$var_info$variable[
        object$srprep$var_info$role == "outcome"
      ]

      # If the outcome isn't included in newdata, use temp placeholder
      if(!any(names(newdata) == outcome_name))
        newdata[[outcome_name]] <- rep(object$srprep$outcome[1], nrow(newdata))

      X <- as.matrix(bake(object$srprep, newdata, everything(), -all_outcomes()))
    } else {
      X <- as.matrix(newdata)
    }
  }



  if(missing(lambda) & length(at) & !length(args$which)) {
    if(at == "cvmin") {
      idx = object$fit$min
    } else {
      idx = which(object$fit$cve < min(object$fit$cve + object$fit$cvse))[1]
    }

    return(predict(object$fit, which = idx, X = X, ...))

  } else if (!missing(lambda)) {

    return(predict(object$fit, lambda = lambda, X=X, ...))

  } else
    predict(object$fit, X=X, ...)
}

#' @rdname predict.sparseR
#' @method coef sparseR
#' @export
coef.sparseR <- function(object, lambda, at = c("cvmin", "cv1se"), ...) {

  at <- match.arg(at)

  args <- list(...)

  if(missing(lambda) & length(at) & !length(args$which)) {
    if(at == "cvmin") {
      idx = object$fit$min
    } else {
      idx = which(object$fit$cve < min(object$fit$cve + object$fit$cvse))[1]
    }

    return(coef(object$fit, which = idx, ...))

  } else if (!missing(lambda)) {
    return(coef(object$fit, lambda = lambda, ...))
  } else
    coef(object$fit, ...)
}


#' Summary of sparseR model coefficients
#' @param object a sparseR object
#' @param lambda a particular value of lambda to predict with
#' @param at a "smart" guess to use for lambda
#' @param ... additional arguments to be passed to summary.ncvreg
#'
#' @method summary sparseR
#'
#' @return an object of class `summary.ncvreg` at specified or smart value of
#'   lambda.
#'
#' @export
summary.sparseR <- function(object, lambda, at = c("cvmin", "cv1se"), ...) {

  at <- match.arg(at)

  args <- list(...)

  if(missing(lambda) & length(at) & !length(args$which)) {
    if(at == "cvmin") {
      idx = object$fit$min
    } else {
      idx = which(object$fit$cve < min(object$fit$cve + object$fit$cvse))[1]
    }

    return(summary(object$fit$fit, which = idx, ...))

  } else if (!missing(lambda)) {
    return(summary(object$fit$fit, lambda = lambda, ...))
  } else
    summary(object$fit$fit, ...)
}
