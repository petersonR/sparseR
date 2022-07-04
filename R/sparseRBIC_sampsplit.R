#' Sample split procedure for stepwise regression
#'
#' Runs multiple on models selection procedures using RBIC to
#' achieve valid inferential results post-selection
#'
#' @param srbic_fit An object fitted by sparseRBIC_step
#' @param S Number of splitting iterations
#' @param quiet Should the display of a progress bar be silenced?
#'
#' @importFrom dplyr full_join
#' @return a list containing:
#'
#' \item{results}{a tibble containing coefficients, p-values, selection pct}
#' \item{splits}{a tibble of different split-based coefficients}
#'
#' @export

sparseRBIC_sampsplit <- function(srbic_fit, S = 100, quiet = FALSE) {

  if(!quiet)
    message("Note: sparseRBIC_sampsplit is currently experimental and may not behave as expected.")

  n <- nrow(srbic_fit$data)

  phi <- list()[1:S]

  pb<-progress_estimated(S)

  for(s in 1:S) {
    idx <- sample(1:n, ceiling(n/2))
    bdata <- bake(srbic_fit$srprep, srbic_fit$data, everything(), -all_outcomes())

    bX <- bdata[idx,]
    by <- srbic_fit$srprep$outcome[idx]

    ## Run stepwise selection with same options as original

    bsrbic_fit <- sparseRBIC_step(
      pre_process = FALSE, model_matrix = bX, y = by,
      direction = "forward",
      family = srbic_fit$info$family,
      cumulative_k = srbic_fit$info$cumulative_k,
      cumulative_poly = srbic_fit$info$cumulative_poly,
      pool = srbic_fit$info$pool,
      poly_prefix = srbic_fit$info$poly_prefix,
      int_sep = srbic_fit$info$int_sep,
      message = FALSE
    )

    ## Fit selected model on new data
    b <- coef(bsrbic_fit$fit)
    active <- gsub("`", "", names(b)[-1])

    bX2 <- bdata[-idx,]
    by2 <- srbic_fit$srprep$outcome[-idx]


    ff <- as.formula(paste0("y~", paste0(active, collapse = " + ")))
    new_fit <- glm(ff, data = data.frame(bX2, y = by2))
    summ <- summary(new_fit)$coef

    phi[[s]] <- tibble(names(coef(bsrbic_fit)), summ[,4])
    names(phi[[s]]) <- c("Coefficient", paste0("s", s))
    if(!quiet)
      pb$tick()$print()
  }

  phi <- Reduce("my_join", phi)

  phi_hat <- as.matrix(phi[,-1]) %>%
    apply(1, median)

  # phi_ci <- as.matrix(phi[,-1]) %>%
  #   apply(1, quantile, probs= c(0.025, .975)) %>%
  #   t() %>%
  #   as.data.frame()

  bootstraps <- t(phi[,-1])
  colnames(bootstraps) <- gsub("`","",phi[[1]])
  bootstraps <- as_tibble(bootstraps, .name_repair = "minimal")

  list(
    results = tibble(
      Coefficient = gsub("`","",phi[[1]]),
      Median_p = phi_hat,
      Mean_p = as.matrix(phi[,-1]) %>%
        apply(1, mean),
      Gmean_p = as.matrix(phi[,-1]) %>%
        apply(1, function(x) exp(mean(log(x)))),
      # Bootstrap_SE = get_bse(Ni, phi, B, phi_hat),
      # p = 2 * pnorm(-abs(Bootstrap_Mean / Bootstrap_SE)),
      "% selected" = 100*apply(phi[,-1], 1, function(x) mean(x != 1))
    ),
    splits = bootstraps
  )
}

my_join <- function(x, y) {
  xy <- dplyr::full_join(x,y, by = "Coefficient")
  xy[is.na(xy)] <- 1
  xy
}

get_bse <- function(Ni, phi, B, phi_hat) {
  ## Calculate Efron smoothed SD
  sd1 <- sqrt(((t(Ni %*% t(as.matrix(phi[,-1]) - phi_hat)) / B) ^ 2) %*% rep(1, nrow(Ni)) )
  as.vector(sd1)
}

