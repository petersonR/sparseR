## Ranked sparsity forward stepwise selection with RBIC
#' Fit a ranked-sparsity model with forward stepwise RBIC (experimental)
#'
#' @param formula Names of the terms
#' @param data Data
#' @param family The family of the model
#' @param k The maximum order of interactions to consider
#' @param poly The maximum order of polynomials to consider
#' @param sequential Should the main effects be considered first, orders
#' sequentially added/considered?
#' @param hier Should hierarchy be enforced (weak or strong)? Must be set
#' with sequential == TRUE (see details)
#' @param cumulative_k Should penalties be increased cumulatively as order
#'   interaction increases?
#' @param cumulative_poly Should penalties be increased cumulatively as order
#'   polynomial increases?
#' @param pool Should interactions of order k and polynomials of order k+1 be
#'   pooled together for calculating the penalty?
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
#' @param trace Should intermediate results of model selection process be output
#' @param ic The information criterion to use
#' @param ia_formula formula to be passed to step_interact via terms argument
#' @param message should experimental message be suppressed
#' @param ... additional arguments for running stepwise selection
#'
#' @details
#'
#' This function mirrors `sparseR` but uses stepwise selection guided by RBIC.
#'
#' Additionally, setting `cumulative_poly` or `cumulative_k` to `TRUE` increases
#' the penalty cumulatively based on the order of either polynomial or
#' interaction.
#'
#' The `hier` hierarchy enforcement will only work if `sequential == TRUE`, and
#' notably will only consider the "first gen" hierarchy, that is, that all
#' main effects which make up an interaction are already in the model. It
#' is therefore possible for a third order interaction (x1:x2:x3) to
#' enter a model without x1:x2 or x2:x3, so long as x1, x2, and x3 are all
#' in the model.
#'
#' The options that can be passed to `pre_proc_opts` are:
#' - knnImpute (should
#' missing data be imputed?)
#' - scale (should data be standardized)?
#' - center
#' (should data be centered to the mean or another value?)
#'  - otherbin (should
#' factors with low prevalence be combined?)
#' - none (should no preprocessing be
#' done? can also specify a null object)
#'
#' The options that can be passed to `extra_opts` are:
#' - centers (named numeric
#' vector which denotes where each covariate should be centered)
#' - center_fn
#' (alternatively, a function can be specified to calculate center such as `min`
#' or `median`)
#' - freq_cut, unique_cut (see ?step_nzv - these get used by the
#' filtering steps)
#' - neighbors (the number of neighbors for knnImpute)
#' - one_hot (see ?step_dummy), this defaults to cell-means coding which can be
#' done in regularized regression (change at your own risk)
#' - raw (should polynomials not be orthogonal? defaults to true because variables are
#' centered and scaled already by this point by default)
#'
#' @md
#'
#' @return an object of class `sparseRBIC` containing the following:
#'
#' \item{fit}{the final fit object}
#' \item{srprep}{a `recipes` object used to prep the data}
#' \item{pen_info}{coefficient-level variable counts, types + names}
#' \item{data}{the (unprocessed) data}
#' \item{family}{the family argument (for non-normal, eg. poisson)}
#' \item{info}{a list containing meta-info about the procedure}
#' \item{stats}{the IC for each fit and respective terms included}
#'
#' @export

sparseRBIC_step <- function(formula, data, family = c("gaussian", "binomial", "poisson"),
                            k = 1, poly = 1, ic = c("RBIC", "RAIC", "BIC", "AIC", "EBIC"),
                            hier = c("strong", "weak", "none"), sequential = (hier[1] != "none"),
                            cumulative_k = FALSE,
                            cumulative_poly = TRUE,
                            pool = FALSE,
                            ia_formula = NULL,
                            pre_process = TRUE, model_matrix = NULL, y = NULL,
                            poly_prefix = "_poly_", int_sep = "\\:",
                            pre_proc_opts = c("knnImpute", "scale", "center", "otherbin", "none"),
                            filter = c("nzv", "zv"), extra_opts = list(),
                            trace = 0, message = TRUE, ...) {

  if(message)
    message("Note: sparseRBIC_step is currently experimental and may not behave as expected.")

  ## Check and validate arguments
  if(is.null(pre_proc_opts))
    pre_proc_opts <- "none"
  pre_proc_opts <- match.arg(pre_proc_opts, several.ok = TRUE)
  filter <- match.arg(filter)
  family <- match.arg(family)
  hier <- match.arg(hier)
  ic <- match.arg(ic)


  if(hier != "none" & !sequential)
    stop("hierarchy can only be enforced if sequential == TRUE")


  if(pre_process) {
    ## Create and prep model matrix function (preprocessing)
    srprep <- sparseR_prep(formula=formula, data=data, k=k, poly=poly, pre_proc_opts=pre_proc_opts,
                           filter=filter, extra_opts=extra_opts, ia_formula = ia_formula)

    ## "Bake" recipe using preprocess steps (this can be done on test data too later)
    X <- bake(srprep, data, everything(), - all_outcomes())

    # Extract outcome
    if(!is.null(srprep$outcome)) {
      y <- srprep$outcome
    } else {
      stop("formula must have outcome for sparseR")
    }

  } else {

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

  if(ic == "RBIC") {
    ic <- RBIC
  } else if(ic == "RAIC"){
    ic <- RAIC
  } else if(ic == "AIC"){
    ic <- myAIC
  } else if(ic == "BIC") {
    ic <- myBIC
  } else if(ic == "EBIC") {
    ic <- EBIC
  }


  ## Get the correct penalty scaling for each covariate
  pen_info <- get_penalties(
    names(X), poly, pool = pool, cumulative_poly = cumulative_poly,
    cumulative_k = cumulative_k, gamma = 1, poly_prefix = poly_prefix,
    int_sep = int_sep
  )

  if((hier == "none" & !sequential) | (k == 0 & poly == 1)) {
    fit <- glm(y ~ 1, data = X, family = family)
    fit_full <- terms(y ~ ., data = X)

    pen_info <- get_penalties(names(X), poly = poly, gamma = 1)

    fit_rbic <- run_steps(object = fit, direction = "forward", scope = list(lower =fit, upper = fit_full),
                          penalty = ic, pen_info = pen_info, varnames = names(X), trace = trace,
                          keep = function(x,...) list(aic = x$aic, terms = attr(x$terms, "term.labels")))

    stats <- fit_rbic$keep
    rownames(stats)[1] <- "IC"

    info <- list(k = k, poly = k, cumulative_k = cumulative_k,
                 cumulative_poly = cumulative_poly, pool = pool,
                 pre_process = pre_process, model_matrix = model_matrix, y = y,
                 poly_prefix = poly_prefix, int_sep = int_sep)

    res <- list(family = family,
                fit = fit_rbic, srprep = srprep, pen_info = pen_info,
                data = data, info = info, stats = stats)

    class(res) <- "sparseRBIC"
    return(res)
  }

  #### Enforcing hierarchy and Sequential
  main_effects_formula <-
    as.formula(paste0("y~", paste0(names(X)[pen_info$vartypes == "Main effect"], collapse = " + ")))
  fit <- glm(y ~ 1, data = X, family = family)
  fit_me <- terms(main_effects_formula, data = X)

  fit_rbic_temp <- run_steps(
    object = fit, direction = "forward", scope = list(lower =fit, upper = fit_me),
    penalty = ic, pen_info = pen_info, varnames = names(X), trace = trace,
    keep = function(x,...) list(aic = x$aic, terms = attr(x$terms, "term.labels"))
  )

  stats <- fit_rbic_temp$keep


  for(order in 1:max(poly-1,k)) {

    matches <-
      grepl(paste("Order", order, "interaction"), pen_info$vartypes) |
      grepl(paste("Order", order+1, "polynomial"), pen_info$vartypes)

    if(!any(matches)) {
      break
    }

    order_k_formula <- terms(as.formula(
        paste0("y~ `", paste0(names(X)[matches], collapse = "` + `"), "`")),
        data = data
      )

    ## To enforce (strong, weak or none) hierarchy for interactions
    vars <- strsplit(gsub("`", "", attributes(order_k_formula)$term.labels), split = int_sep)
    keep_idx <- rep(TRUE, length(vars))

    int_idx <- sapply(vars, length) > 1

    if(k > 0) {
      keep_idx[int_idx] <- sapply(vars[int_idx], function(v) {
        vnames <- gsub(paste0(poly_prefix, "1"), "",names(fit_rbic_temp$coef))
        if(hier == "weak")
          return(any(vnames %in% v))
        if(hier == "strong")
          return(all(v %in% vnames))
        TRUE
      })
    }

    ## To enforce (strong, weak or none) hierarchy for polynomials
    ## (looks for polynomial of order k)
    if(poly > 1 & any(!int_idx)) {
      keep_idx[!int_idx] <- sapply(vars[!int_idx], function(v) {
        vname <- paste0(substring(v, 0, nchar(v) - 1), order)
        if(hier == "weak" | hier == "strong")
          return(any(names(fit_rbic_temp$coef) == vname))
        TRUE
      })
    }

    if(any(unlist(keep_idx))) {
      allowed_terms <- as.formula(
        paste0("y~.+", paste0(attributes(order_k_formula)$term.labels[keep_idx], collapse = " + ")))

    } else
      allowed_terms <- y~.

    updated_formula <-  update(terms(fit_rbic_temp), allowed_terms)

    fit_rbic_temp <- run_steps(object = fit_rbic_temp, direction = "forward",
                             scope = list(lower =fit_rbic_temp, upper = updated_formula),
                             penalty = ic, pen_info = pen_info, varnames = names(X),
                             trace = trace,
                             keep = function(x,...) list(aic = x$aic, terms = attr(x$terms, "term.labels")))
    stats <- cbind(stats, fit_rbic_temp$keep[,-1])
  }

  rownames(stats)[1] <- "IC"
  fit_rbic <- fit_rbic_temp

  info <- list(k = k, poly = k, sequential = sequential,
               hier = hier,
               cumulative_k = cumulative_k,
               cumulative_poly = cumulative_poly, pool = pool,
               pre_process = pre_process, model_matrix = model_matrix, y = y,
               poly_prefix = poly_prefix, int_sep = int_sep)

  res <- list(fit = fit_rbic, srprep = srprep, pen_info = pen_info,
              data = data, family = family,
              info = info, stats = t(stats))

  class(res) <- "sparseRBIC"
  res
}

#' @method print sparseRBIC
#' @export
print.sparseRBIC <- function(x, ...) {
  print(x$fit)
}

#' @method coef sparseRBIC
#' @export

coef.sparseRBIC <- function(object, ...) {
  coef(object$fit)
}

#' @method predict sparseRBIC
#' @export

predict.sparseRBIC <- function(object, newdata, ...) {

  ## If newdata is not specified...
  if(missing(newdata)) {
    if(!is.null(object$srprep)) {
      X <- bake(object$srprep, object$data, everything(), -all_outcomes())
    } else {
      X <- NULL
    }
  } else {
    if(!is.null(object$srprep)) {
      X <- bake(object$srprep, newdata, everything(), -all_outcomes())
    } else {
      X <- newdata
    }
  }
  predict(object$fit, newdata = X, ...)
}

#' @method summary sparseRBIC
#' @export

summary.sparseRBIC <- function(object, ...) {
  message("P-values have **not** been corrected for multiple comparisons\n",
          "Consider sparseRBIC_sampsplit() or sparseRBIC_bootstrap()")
  summary(object$fit, ...)
}
