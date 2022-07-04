#' Centering numeric data to a value besides their mean
#'
#' `step_center_to` generalizes `step_center` to allow for a different function
#' than the `mean` function to calculate centers. It creates a *specification*
#' of a recipe step that will normalize numeric data to have a 'center' of zero.
#'
#' @param recipe A recipe object. The step will be added to the sequence of
#'   operations for this recipe.
#' @param ... One or more selector functions to choose which variables are
#'   affected by the step. See [selections()] for more details. For the `tidy`
#'   method, these are not currently used.
#' @param role Not used by this step since no new variables are created.
#' @param trained A logical to indicate if the quantities for preprocessing have
#'   been estimated.
#' @param centers A named numeric vector of centers. This is `NULL` until
#'   computed by [prep.recipe()] (or it can be specified as a named
#'   numeric vector as well?).
#' @param na_rm A logical value indicating whether `NA` values should be removed
#'   during computations.
#' @param skip A logical. Should the step be skipped when the recipe is baked by
#'   [bake.recipe()]? While all operations are baked when [prep.recipe()] is
#'   run, some operations may not be able to be conducted on new data (e.g.
#'   processing the outcome variable(s)). Care should be taken when using `skip
#'   = TRUE` as it may affect the computations for subsequent operations
#' @param id A character string that is unique to this step to identify it.
#' @param center_fn a function to be used to calculate where the center should be
#'
#' @return An updated version of `recipe` with the new step added to the
#'   sequence of existing steps (if any). For the `tidy` method, a tibble with
#'   columns `terms` (the selectors or variables selected) and `value` (the
#'   centers).
#'
#' @import recipes
#' @importFrom dplyr as_tibble
#'
#' @keywords datagen
#' @concept preprocessing
#' @concept normalization_methods
#' @export
#' @details Centering data means that the average of a variable is subtracted
#'   from the data. `step_center_to` estimates the variable centers from the
#'   data used in the `training` argument of `prep.recipe`. `bake.recipe` then
#'   applies the centering to new data sets using these centers.
#'
#' @examples
#' data(biomass, package = "modeldata")
#'
#' biomass_tr <- biomass[biomass$dataset == "Training",]
#' biomass_te <- biomass[biomass$dataset == "Testing",]
#'
#' rec <- recipes::recipe(
#'  HHV ~ carbon + hydrogen + oxygen + nitrogen + sulfur,
#'  data = biomass_tr)
#'
#' center_trans <- rec %>%
#'   step_center_to(carbon, contains("gen"), -hydrogen)
#'
#' center_obj <- recipes::prep(center_trans, training = biomass_tr)
#'
#' transformed_te <- recipes::bake(center_obj, biomass_te)
#'
#' biomass_te[1:10, names(transformed_te)]
#' transformed_te
#'
#' recipes::tidy(center_trans)
#' recipes::tidy(center_obj)
#' @seealso [recipe()] [prep.recipe()] [bake.recipe()]
step_center_to <-
  function(recipe,
           ...,
           role = NA,
           trained = FALSE,
           centers = NULL,
           center_fn = mean,
           na_rm = TRUE,
           skip = FALSE,
           id = rand_id("center_to")) {
    add_step(
      recipe,
      step_center_to_new(
        terms = ellipse_check(...),
        trained = trained,
        role = role,
        centers = centers,
        center_fn = center_fn,
        na_rm = na_rm,
        skip = skip,
        id = id
      )
    )
  }

## Initializes a new object
step_center_to_new <-
  function(terms, role, trained, centers, center_fn, na_rm, skip, id) {
    step(
      subclass = "center_to",
      terms = terms,
      role = role,
      trained = trained,
      centers = centers,
      center_fn = center_fn,
      na_rm = na_rm,
      skip = skip,
      id = id
    )
  }

#' @method prep step_center_to
#' @export
#'
prep.step_center_to <- function(x, training, info = NULL, ...) {
  col_names <- recipes_eval_select(x$terms, training, info = info)
  check_type(training[, col_names])

  if(is.null(x$centers)) {
    centers <-
      vapply(training[, col_names], x$center_fn, c(center = 0), na.rm = x$na_rm)
    attr(centers, "custom") <- FALSE
  } else
    centers <- x$centers

  step_center_to_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    centers = centers,
    center_fn = x$center_fn,
    na_rm = x$na_rm,
    skip = x$skip,
    id = x$id
  )
}

#' @method bake step_center_to
#' @export
bake.step_center_to <- function(object, new_data, ...) {
  res <-
    sweep(as.matrix(new_data[, names(object$centers)]), 2, object$centers, "-")
  if (is.matrix(res) && ncol(res) == 1)
    res <- res[, 1]
  new_data[, names(object$centers)] <- res
  as_tibble(new_data)
}

#' @method print step_center_to
#' @export
print.step_center_to <-
  function(x, width = max(20, options()$width - 30), ...) {

    custom <- length(attr(x$centers, "custom")) && attr(x$centers, "custom") |
      length(attr(x$center_fn, "custom")) && attr(x$center_fn, "custom")


    if(custom)
      cat("Centering to custom value for ",  sep = "")
    else
      cat("Centering to mean for ", sep = "")
    printer(names(x$centers), x$terms, x$trained, width = width)
    invisible(x)
  }


#' @rdname step_center_to
#' @param x A `step_center_to` object.
#' @method tidy step_center_to
#' @importFrom rlang na_dbl
#' @export
tidy.step_center_to <- function(x, ...) {
  if (is_trained(x)) {
    res <- tibble(terms = names(x$centers),
                  value = x$centers)
  } else {
    term_names <- sel2char(x$terms)
    res <- tibble(terms = term_names,
                  value = na_dbl)
  }
  res$id <- x$id
  res
}
