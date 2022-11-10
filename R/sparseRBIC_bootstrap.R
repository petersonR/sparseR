#' Bootstrap procedure for stepwise regression
#'
#' Runs bootstrap on models selection procedure using RBIC to
#' find bootstrapped standard error (smoothed, see Efron 2014)
#' as well as selection percentage across candidate variables.
#' (experimental)
#'
#' @param srbic_fit An object fitted by sparseRBIC_step
#' @param B Number of bootstrap samples
#' @param quiet Should the display of a progress bar be silenced?
#' @importFrom dplyr progress_estimated tibble
#'
#' @return a list containing:
#'
#' \item{results}{a tibble containing coefficients, p-values, selection pct}
#' \item{bootstraps}{a tibble of bootstrapped coefficients}
#'
#' @export

sparseRBIC_bootstrap <- function(srbic_fit, B = 100, quiet = FALSE) {

 if(!quiet)
   message("Note: sparseRBIC_bootstrap is currently experimental and may not behave as expected.")


  ## Need to keep track of how many times each observation shows up
  Ni <- matrix(0,nrow = nrow(srbic_fit$data), ncol = B)
  mode(Ni) <- "integer"

  phi <- list()[1:B]

  if(!quiet)
    pb<-dplyr::progress_estimated(B)

  for(b in 1:B) {
    idx <- sample(1:nrow(srbic_fit$data), replace = TRUE)
    t <- table(idx)

    Ni[as.numeric(names(as.array(t))), b] <- as.numeric(t)

    bdata <- srbic_fit$data[idx,]

    bX <- bake(srbic_fit$srprep, bdata, everything(), -all_outcomes())

    # Extract outcome
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

    s <- summary(bsrbic_fit$fit)$coef

    phi[[b]] <- tibble(names(coef(bsrbic_fit)), s[,4])
    names(phi[[b]]) <- c("Coefficient", paste0("b", b))
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
      Mean_p = as.matrix(phi[,-1]) %>%
        apply(1, mean),
      Gmean_p = as.matrix(phi[,-1]) %>%
        apply(1, function(x) exp(mean(log(x)))),
      # Bootstrap_SE = get_bse(Ni, phi, B, phi_hat),
      # p = 2 * pnorm(-abs(Bootstrap_Mean / Bootstrap_SE)),
      "% selected" = 100*apply(phi[,-1], 1, function(x) mean(x != 1))
    ),
    bootstraps = bootstraps
  )
}

get_bse <- function(Ni, phi, B, phi_hat) {
  ## Calculate Efron smoothed SD
  sd1 <- sqrt(((t(Ni %*% t(as.matrix(phi[,-1]) - phi_hat)) / B) ^ 2) %*% rep(1, nrow(Ni)) )
  as.vector(sd1)
}

