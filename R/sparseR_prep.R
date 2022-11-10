#' Preprocess & create a model matrix with interactions + polynomials
#' @param formula A formula of the main effects + outcome of the model
#' @param data A required data frame or tibble containing the variables in
#'   \code{formula}
#' @param k Maximum order of interactions to *numeric* variables
#' @param pre_proc_opts A character vector specifying methods for preprocessing
#'   (see details)
#' @param ia_formula formula to be passed to step_interact (for interactions,
#'   see details)
#' @param filter which methods should be used to filter out variables with
#'   (near) zero variance? (see details)
#' @param poly the maximum order of polynomials to consider
#' @param extra_opts extra options to be used for preprocessing
#' @param family family passed from sparseR
#' @md
#' @details
#'
#' The pre_proc_opts acts as a wrapper for the corresponding procedures in the
#' \code{recipes} package. The currently supported options that can be passed to
#' pre_proc_opts are: knnImpute: Should k-nearest-neighbors be performed (if
#' necessary?) scale: Should variables be scaled prior to creating interactions
#' (does not scale factor variables or dummy variables) center: Should variables
#' be centered (will not center factor variables or dummy variables ) otherbin:
#'
#' \code{ia_formula} will by default interact all variables with each other up
#' to order k. If specified, ia_formula will be passed as the `terms` argument
#' to \code{recipes::step_interact}, so the help documentation for that function
#' can be investigated for further assistance in specifying specific
#' interactions.
#'
#' The methods specified in filter are important; filtering is necessary to cut
#' down on extraneous polynomials and interactions (in cases where they really
#' don't make sense). This is true, for instance, when using dummy variables in
#' polynomials , or when using interactions of dummy variables that relate to
#' the same categorical variable.
#'
#' @return an object of class `recipe`; see [recipes::recipe()]
#'
#' @importFrom dplyr all_of
#' @importFrom recipes recipe step_nzv step_zv has_role all_predictors all_outcomes
#'   prep update_role step_interact step_naomit all_numeric step_scale step_other
#'   step_impute_knn step_naomit step_rm step_dummy step_poly
#'
#' @export
sparseR_prep <- function(formula, data, k = 1, poly = 1,
                         pre_proc_opts = c("knnImpute", "scale", "center", "otherbin", "none"),
                         ia_formula = NULL,
                         filter = c("nzv", "zv"), extra_opts = list(),
                         family = "gaussian") {


  stopifnot((as.integer(k) == k) && k >= 0)
  stopifnot((as.integer(poly) == poly) && poly > 0)

  # Create recipe
  rec_obj <- recipe(formula, data = data)

  # filter out survival outcomes
  outcome <- NULL
  if(family == "coxph" ) {
    if(any(rec_obj$var_info$role == "outcome")) {
      outcome <- bake(prep(rec_obj, data), all_outcomes(), new_data =data)[[1]]
      rec_obj <- rec_obj %>% step_rm(all_outcomes())
    }
  }

  # must filter out censored predictors as well or will crash
  rec_obj <- rec_obj %>% step_rm(has_type("censored"), has_type("surv"))

  ## Filter out near-zero variance main effects
  if("nzv" %in% filter) {
    if(!length(extra_opts$unique_cut))
      extra_opts$unique_cut <- 10
    if(!length(extra_opts$freq_cut))
      extra_opts$freq_cut <- 95/5

    rec_obj <- rec_obj %>%
      step_nzv(all_predictors(), unique_cut = extra_opts$unique_cut,
               freq_cut = extra_opts$freq_cut)
  } else if("zv" %in% filter) {
    rec_obj <- rec_obj %>% step_zv(all_predictors())
  }

  # Label any dummy variables correctly (requires prepping)
  p_early <- prep(rec_obj, data)
  dummies <- p_early$template %>%
    sapply(function(x) all(x[!is.na(x)] %in% c(0,1)))
  dummies <- names(dummies)[dummies]
  dummies <- dummies[dummies != rec_obj$var_info$variable[rec_obj$var_info$role == "outcome"]]

  if(length(dummies)) {
    rec_obj <- rec_obj %>%
      update_role(all_of(dummies), new_role = "dummy")
  }

  # Center
  p_early <- prep(rec_obj, data)
  has_numeric_predictor <- any(
    vapply(
      X = p_early$term_info$type[p_early$term_info$role == "predictor"],
      FUN = function(x) "numeric" %in% x,
      FUN.VALUE = logical(1)
    )
  )

  if("center" %in% pre_proc_opts & has_numeric_predictor) {
    if(!length(extra_opts$centers)) {
      extra_opts$centers <- NULL
    } else
      attr(extra_opts$centers, "custom") <- TRUE

    if(!length(extra_opts$center_fn)) {
      extra_opts$center_fn <- mean
      attr(extra_opts$center_fn, "custom") <- FALSE
    } else
      attr(extra_opts$center_fn, "custom") <- TRUE

    rec_obj <- rec_obj %>%
      step_center_to(all_numeric(), -has_role("dummy"), - all_outcomes(),
                     centers = extra_opts$centers,
                     center_fn = extra_opts$center_fn)
  }

  # Scale
  if("scale" %in% pre_proc_opts & has_numeric_predictor) {
    rec_obj <- rec_obj %>%
      step_scale(all_numeric(), -has_role("dummy"), -all_outcomes())
  }

  # Bin "other" factors if necessary
  if("otherbin" %in% pre_proc_opts)
    rec_obj <- rec_obj %>% step_other(all_nominal(), -all_outcomes())

  # KNN imputation
  if("knnImpute" %in% pre_proc_opts) {
    if(!length(extra_opts$neighbors))
      extra_opts$neighbors <- 5
    rec_obj <- rec_obj %>%
      step_impute_knn(all_predictors(), has_role("dummy"),
                     neighbors = extra_opts$neighbors,
                     impute_with = imp_vars(all_predictors(), has_role("dummy")))
  }


  # If KNN failed, use step_na
  rec_obj <- rec_obj %>%
    step_naomit(everything(), skip = FALSE)

  # Create dummy variables if necessary
  if(!length(extra_opts$one_hot))
    extra_opts$one_hot <- TRUE


  if(family != "coxph" & any(rec_obj$var_info$role == "outcome")) {
    outcome <- bake(prep(rec_obj, data), all_outcomes(), new_data = data)[[1]]
    rec_obj <- rec_obj %>% step_rm(all_outcomes())
  }


  rec_obj_early <- rec_obj %>%
    step_dummy(all_nominal(),
               one_hot = extra_opts$one_hot, role = "dummy") %>%
    prep(data)

  rec_obj <- rec_obj_early

  # Create interactions if necessary
  prec_obj <- prep(rec_obj, data = data)
  has_predictors <- any(prec_obj$term_info$role == "predictor")

  if(k > 0) {

    has_dummies <- any(prec_obj$term_info$role == "dummy")

    if(is.null(ia_formula)) {
      if(has_dummies & has_predictors) {
        form_var <- "(all_predictors() + has_role('dummy'))"
      } else if(has_dummies) {
        form_var <- "has_role('dummy')"
      } else if(has_predictors) {
        form_var <- "all_predictors()"
      } else {
        stop("No predictors left to make interactions (check if filtered?)")
      }

      ## Need to specify formula differently if dummy vars are included
      ia_formula <- as.formula(paste("~",
                                     paste(rep(form_var, k + 1), collapse = " * "))
      )
    } else if (!("formula" %in% class(ia_formula)))
      stop("ia_formula must be of class formula, see ?recipes::step_interact")

    rec_obj  <- rec_obj %>%
      step_interact(terms = ia_formula, sep = ":", role = "interaction")
  } else if (!is.null(ia_formula)) {
    stop("k must be > 0 for specific ia_formula")
  }

  if(poly > 1 & has_predictors) {
    if(!length(extra_opts$raw))
      extra_opts$raw <- TRUE

    rec_obj <- rec_obj %>%
      step_poly(all_predictors(), degree = poly,
                options = list(raw = extra_opts$raw))
  }

  ## Filter out near-zero variance interaction/polynomial effects
  if("nzv" %in% filter) {
    rec_obj <- rec_obj %>%
      step_nzv(everything(), unique_cut = extra_opts$unique_cut,
               freq_cut = extra_opts$freq_cut)
  } else if("zv" %in% filter) {
    rec_obj <- rec_obj %>% step_zv(everything())
  } else
    warning("Recommend at least filtering zv for extraneous interaction effects")

  ## Prep recipe
  final_rec <- rec_obj %>%
    prep(data)

  ## Fix tr_info for earlier prepped version
  final_rec$tr_info <- rec_obj_early$tr_info

  # Add meta-information
  final_rec$pre_proc_opts <- pre_proc_opts
  final_rec$filter <- filter
  final_rec$extra_opts <- extra_opts
  final_rec$outcome <- outcome

  final_rec
}

