#' Plot relevant properties of sparseR objects
#'
#' @param x a `sparseR` object
#' @param plot_type should the solution path, CV results, or both be plotted?
#' @param cols option to specify color of groups
#' @param log.l should the x-axis (lambda) be logged?
#' @param ... extra plotting options
#'
#' @method plot sparseR
#'
#' @return nothing returned
#'
#' @export
plot.sparseR <- function(x, plot_type = c("both", "cv", "path"), cols = NULL, log.l = TRUE, ...) {

  plot_type <- match.arg(plot_type)

  vline <- c(
    x$fit$lambda[which.min(x$fit$cve)],
    x$fit$lambda[which(x$fit$cve < min(x$fit$cve + x$fit$cvse))[1]]
  )

  if(log.l)
    vline <- log(vline)


  if("both" %in% plot_type | "cv" %in% plot_type) {
    plot(x$fit,  log.l=log.l, ...)
    abline(v = vline, lty = c(2), col =c("slategrey", "orange3"))
  }

  if("both" %in% plot_type | "path" %in% plot_type) {

    if(!length(cols)) {

      ncols <- length(unique(x$results$Vartype))
      cols <- factor(x$results$Vartype, labels = pal(ncols))

    } else
      cols <- factor(x$results$Vartype, levels = names(cols), labels = cols)


    ## Need to get index for coloring (to match plot.ncvreg)
    YY <- if (length(x$fit$fit$penalty.factor) == nrow(x$fit$fit$beta)) {
      coef(x$fit$fit)
    } else coef(x$fit$fit)[-1, , drop = FALSE]


    # only include those coefficients which are penalized
    ind <- which(x$fit$fit$penalty.factor != 0)

    plot(x$fit$fit, log.l = log.l, col = as.character(cols[ind]))
    abline(v = vline, lty = c(2), col =c("slategrey", "orange3"))
    legend("topleft", legend = levels(factor(x$results$Vartype)),
           bty = "n", col = levels(cols), lwd = 2)

  }
  return(invisible(NULL))
}

## Adapted from Patrick Breheny's ncvreg code
pal <- function (n, alpha = 1) {
  if (n == 2) {
    val <- hcl(seq(15, 375, len = 4), l = 60, c = 150, alpha = alpha)[c(1,
                                                                        3)]
  }
  else val <- hcl(seq(15, 375, len = n + 1), l = 60, c = 150,
                  alpha = alpha)[1:n]
  val
}

#' Plot relevant effects of a sparseR object
#'
#' @rdname effect_plot
#'
#'
#'
#' @param fit a `sparseR` object
#' @param coef_name The name of the coefficient to plot along the x-axis
#' @param at value of lambda to use
#' @param by the variable(s) involved in the (possible) interaction
#' @param by_levels values to cut continuous by variable (defaults to 3 quantiles)
#' @param nn number of points to plot along prediction line
#' @param plot.args list of arguments passed to the plot itself
#' @param resids should residuals be plotted or not?
#' @param legend_location location for legend passed to `legend`
#' @param ... additional arguments
#' @return nothing returned
#' @export

effect_plot <- function(fit, ...) {
  UseMethod("effect_plot")
}

#' @rdname effect_plot
#' @export
#' @method effect_plot sparseR
effect_plot.sparseR <- function(fit, coef_name, at = c("cvmin", "cv1se"),
                        by = NULL, by_levels, nn = 101,
                        plot.args = list(), resids = TRUE,
                        legend_location = "bottomright",
                        ...) {


  if(is.null(fit$srprep))
    stop("effect_plot currently only supports sparseR models specified w. formula")

  outcome_name <- fit$srprep$var_info$variable[fit$srprep$var_info$role == "outcome"]
  yy <- fit$srprep$outcome
  if(inherits(yy, what = "Surv"))
    resids <- FALSE
  pclass <- sapply(fit$data, class)

  ## Send warning if not getting the "full picture" of interaction
  b <- coef(fit, at  = at)
  active <- names(b[b!= 0])
  active_int <- active[grepl("\\:", active)]

  if(any(relevant <- grepl(coef_name, active_int))) {

    if(!is.null(by)) {
      unaccounted <- active_int[relevant & (!grepl(by, active_int))]
      if(length(unaccounted)) {
        warning(paste0(coef_name,
                       " has the following interactions which are not displayed: ",
                       paste0(unaccounted, collapse = ", ")))

      }
    } else {
       warning(paste0(coef_name,
                     " has the following interactions which are not displayed: ",
                     paste0(active_int[relevant], collapse = ", ")))
    }
  }

  ## Set up model matrix with medians for all numeric
  new_data <- fit$data[rep(1, nn),]
  new_data[,pclass == "numeric"] <-
    sapply(fit$data[pclass == "numeric"], function(x) rep(median(x, na.rm = TRUE), nn))

  ## And most prevalent group for all factors
  new_data[,pclass == "factor"] <-
    lapply(fit$data[pclass == "factor"],
           function(x) rep(levels(x)[which.max(table(x))], nn))

  ## If the coef_name is a factor variable, get relevant info
  if(cfact <- any(class(new_data[[coef_name]]) %in% c("factor","character"))) {
    levs <- levels(fit$data[[coef_name]])
    new_data[[coef_name]] <- rep(levs, nn)[1:nn]
    new_data <- new_data[1:length(levs),]
  } else
    new_data[[coef_name]] <- seq(min(fit$data[[coef_name]]), max(fit$data[[coef_name]]), length = nn)

  if(!length(legend_location))
    legend_location <- "bottomright"

  # Validate plot.args
  if(!length(plot.args$ylab))
    plot.args$ylab <- outcome_name

  if(!length(plot.args$xlab))
    plot.args$xlab <- coef_name

  if(!length(plot.args$points.col))
    plot.args$points.col <- "slategrey"

  ### No by variable
  if(!length(by)) {

    preds <- predict(fit, newdata=new_data, at = at)

    ### Numeric, no by variable
    if(!cfact) {
      resid <- yy - predict(fit, at = at) +
        approx(new_data[[coef_name]], preds, xout = fit$data[[coef_name]])$y

      if(!length(plot.args$ylim))
        plot.args$ylim <- range(preds, resid, na.rm = TRUE)

      if(!length(plot.args$xlim))
        plot.args$xlim <- range(new_data[[coef_name]])

      plot(new_data[[coef_name]], preds, type = "l", lwd = 2,
           xlab = plot.args$xlab, ylab = plot.args$ylab,
           ylim = plot.args$ylim, xlim = plot.args$xlim, ...)

      if(resids)
        points(fit$data[[coef_name]], resid, pch = 20, col = plot.args$points.col)

      ### Factor covariate, no by variable
    } else {

      xx <- sapply(1:length(levs), function(x) seq(.75, 1.25, length = nn/length(levs)) + x-1)
      pp <- rep(preds, each = nrow(xx))

      resid <- yy - predict(fit, at = at) +
        as.numeric(as.character(factor(fit$data[[coef_name]], labels = preds)))

      if(!length(plot.args$ylim))
        plot.args$ylim <- range(preds, resid, na.rm = TRUE)

      if(!length(plot.args$xlim))
        plot.args$xlim <- range(xx)

      plot(c(1), c(1), type = "n",
           xlab = plot.args$xlab, ylab = plot.args$ylab,
           ylim = plot.args$ylim, xlim = plot.args$xlim, xaxt = "n", ...)

      xx <- sapply(1:length(levs), function(i) {
        xx <- seq(.7, 1.3, length = nn/length(levs)) + i-1
        lines(xx, rep(preds[i], each= length(xx)), lwd = 2)
        xx
      })

      axis(1, at = 1:length(levs), levs)

      if(resids)
        points(jitter(as.numeric(fit$data[[coef_name]])), resid, pch = 20, col = "slategrey")
    }

    #### IF there IS a by variable
  } else {

    ### Factor by variable
    if(inherits(fit$data[[by]], what = "factor")) {
      by_levels <- levels(fit$data[[by]])
      pt_cols <- as.character(factor(fit$data[[by]], labels = pal(length(by_levels))))

      #### Numeric covariate ####
      if(!cfact) {

        resid <- preds <- list()
        for(b in 1:length(by_levels)) {
          new_data[[by]] <- rep(by_levels[[b]], nn)

          preds[[b]] <- predict(fit, newdata = new_data, at = at)
          if(resids)
            resid[[b]] <- (yy[fit$data[[by]] == by_levels[b]]) -
              predict(fit, at = at, newdata = fit$data[fit$data[[by]] == by_levels[b],]) +
              approx(new_data[[coef_name]], preds[[b]], xout = fit$data[[coef_name]][fit$data[[by]] == by_levels[b]])$y
        }


        # Validate plot.args
        if(!length(plot.args$ylab))
          plot.args$ylab <- outcome_name

        if(!length(plot.args$xlab))
          plot.args$xlab <- coef_name

        if(!length(plot.args$ylim))
          plot.args$ylim <- range(preds, resid)

        if(!length(plot.args$xlim))
          plot.args$xlim <- range(new_data[[coef_name]])

        cols <- pal(length(by_levels))
        plot(new_data[[coef_name]], preds[[1]], type = "n", lwd = 2,
             xlab = plot.args$xlab, ylab = plot.args$ylab,
             ylim = plot.args$ylim, xlim = plot.args$xlim, ...)



        lty <- rep(1, length(by_levels))
        for(i in 1:length(by_levels)) {
          ## Check for overlapping lines
          if(i > 2) {
            overlap <- which(sapply(1:(i-1), function(j) isTRUE(all.equal(preds[[j]], preds[[i]]))))
          } else
            overlap <- NULL

          if(length(overlap))
            lty[i] = 2

          lines(new_data[[coef_name]], preds[[i]], lwd = 2, col = cols[i], lty = lty[i])
        }

        points(fit$data[[coef_name]], unlist(resid), pch = 20, col = pt_cols)


        leg_lab <- paste0(by, ": ", unname(by_levels))
        legend(legend_location, legend = leg_lab, col = cols, lwd = 2, bty = "n", lty = lty)

        ## Factor main covariate
      } else
        stop("Cannot yet handle factor effect for both main covariate and 'by' covariate")

    ## Numeric by variable
    } else if(inherits(fit$data[[by]], what = "numeric")) {


      ## IF it's a numeric variable
      if(!cfact) {
        if(missing(by_levels)) {
          by_levels <- quantile(fit$data[[by]], c(.1, .5, .9), na.rm = TRUE)
        }
        closest_level <- apply(sapply(by_levels, function(x) (abs(x - fit$data[[by]]))), 1, which.min)

        resid <- preds <- list()
        for(b in 1:length(by_levels)) {
          new_data[[by]] <- rep(by_levels[[b]], nn)

          preds[[b]] <- predict(fit, newdata=new_data, at = at)
          resid[[b]] <- yy[closest_level == b] -
            predict(fit, at = at, newdata = fit$data[closest_level == b,]) +
            approx(new_data[[coef_name]], preds[[b]], xout = fit$data[[coef_name]][closest_level == b])$y
        }


        # Validate plot.args
        if(!length(plot.args$ylab))
          plot.args$ylab <- outcome_name

        if(!length(plot.args$xlab))
          plot.args$xlab <- coef_name

        if(!length(plot.args$ylim))
          plot.args$ylim <- range(preds)

        if(!length(plot.args$xlim))
          plot.args$xlim <- range(new_data[[coef_name]])

        cols <- pal(length(by_levels))
        plot(new_data[[coef_name]], preds[[1]], type = "n", lwd = 2,
             xlab = plot.args$xlab, ylab = plot.args$ylab,
             ylim = plot.args$ylim, xlim = plot.args$xlim, ...)

        ll <- lapply(1:length(by_levels), function(i) {
          lines(new_data[[coef_name]], preds[[i]], lwd = 2, col = cols[i])
        })
        if (resids) {
          pp <- lapply(1:length(by_levels), function(i) {
            points(fit$data[[coef_name]][closest_level == i], resid[[i]], pch = 20, col = cols[i])
          })
        }

        leg_lab <- paste0(by, ": ", round(unname(by_levels),2))
        legend(legend_location, legend = leg_lab, col = cols, lwd = 2, bty = "n")

        ## factor covariate by numeric one
      } else {
        stop("Cannot yet handle factor covariate by numeric covariate")
      }
    }
  }
  return(invisible(NULL))
}


#' @rdname effect_plot
#'
#' @export
#' @method effect_plot sparseRBIC
#' @return Nothing (invisible) returned
effect_plot.sparseRBIC <- function(fit, coef_name,
                                by = NULL, by_levels, nn = 101,
                                plot.args = list(), resids = TRUE,
                                legend_location =  "bottomright", ...) {


  if(is.null(fit$srprep))
    stop("effect_plot currently only supports sparseR models specified w. formula")

  outcome_name <- fit$srprep$var_info$variable[fit$srprep$var_info$role == "outcome"]
  pclass <- sapply(fit$data, class)
  yy <- fit$srprep$outcome

  ## Send warning if not getting the "full picture" of interaction
  b <- coef(fit)
  active <- names(b[b!= 0])
  active_int <- active[grepl("\\:", active)]

  if(any(relevant <- grepl(coef_name, active_int))) {

    if(!is.null(by)) {
      unaccounted <- active_int[relevant & (!grepl(by, active_int))]
      if(length(unaccounted)) {
        warning(paste0(coef_name,
                       " has the following interactions which are not displayed: ",
                       paste0(unaccounted, collapse = ", ")))

      }
    } else {
      warning(paste0(coef_name,
                     " has the following interactions which are not displayed: ",
                     paste0(active_int[relevant], collapse = ", ")))
    }
  }

  ## Set up model matrix with medians for all numeric
  new_data <- fit$data[rep(1, nn),]
  new_data[pclass == "numeric"] <-
    sapply(fit$data[pclass == "numeric"], function(x) rep(median(x, na.rm = TRUE), nn))

  ## And most prevalent group for all factors
  new_data[pclass == "factor"] <-
    sapply(fit$data[pclass == "factor"], function(x) rep(levels(x)[which.max(table(x))], nn))

  ## If the coef_name is a factor variable, get relevant info
  if(cfact <- any(class(fit$data[[coef_name]]) %in% c("factor","character"))) {
    levs <- levels(fit$data[[coef_name]])
    new_data[[coef_name]] <- rep(levs, nn)[1:nn]
    new_data <- new_data[1:length(levs),]
  } else
    new_data[[coef_name]] <- seq(min(fit$data[[coef_name]], na.rm = TRUE),
                                 max(fit$data[[coef_name]], na.rm = TRUE), length = nn)

  # Validate plot.args
  if(!length(plot.args$ylab))
    plot.args$ylab <- outcome_name

  if(!length(plot.args$xlab))
    plot.args$xlab <- coef_name

  if(!length(plot.args$points.col))
    plot.args$points.col <- "slategrey"

  ### No by variable
  if(!length(by)) {

    preds <- predict(fit, newdata=new_data)

    ### Numeric, no by variable
    if(!cfact) {
      resid <- yy - predict(fit) +
        approx(new_data[[coef_name]], preds, xout = fit$data[[coef_name]])$y

      if(!length(plot.args$ylim))
        plot.args$ylim <- range(preds, resid, na.rm = TRUE)

      if(!length(plot.args$xlim))
        plot.args$xlim <- range(new_data[[coef_name]])

      plot(new_data[[coef_name]], preds, type = "l", lwd = 2,
           xlab = plot.args$xlab, ylab = plot.args$ylab,
           ylim = plot.args$ylim, xlim = plot.args$xlim, ...)

      if(resids)
        points(fit$data[[coef_name]], resid, pch = 20, col = plot.args$points.col)

      ### Factor covariate, no by variable
    } else {

      xx <- sapply(1:length(levs), function(x) seq(.75, 1.25, length = nn/length(levs)) + x-1)
      pp <- rep(preds, each = nrow(xx))

      resid <- yy - predict(fit) +
        as.numeric(as.character(factor(fit$data[[coef_name]], labels = preds)))

      if(!length(plot.args$ylim))
        plot.args$ylim <- range(preds, resid, na.rm=TRUE)

      if(!length(plot.args$xlim))
        plot.args$xlim <- range(xx)

      plot(c(1), c(1), type = "n",
           xlab = plot.args$xlab, ylab = plot.args$ylab,
           ylim = plot.args$ylim, xlim = plot.args$xlim, xaxt = "n", ...)

      xx <- sapply(1:length(levs), function(i) {
        xx <- seq(.7, 1.3, length = nn/length(levs)) + i-1
        lines(xx, rep(preds[i], each= length(xx)), lwd = 2)
        xx
      })

      axis(1, at = 1:length(levs), levs)

      if(resids)
        points(jitter(as.numeric(fit$data[[coef_name]])), resid, pch = 20, col = "slategrey")
    }

    #### IF there IS a by variable
  } else {

    ### Factor by variable
    if(inherits(fit$data[[by]], what = "factor")) {
      by_levels <- levels(fit$data[[by]])
      pt_cols <- as.character(factor(fit$data[[by]], labels = pal(length(by_levels))))

      #### Numeric covariate ####
      if(!cfact) {

        resid <- preds <- list()
        for(b in 1:length(by_levels)) {
          new_data[[by]] <- rep(by_levels[[b]], nn)

          preds[[b]] <- predict(fit, newdata = new_data)
          resid[[b]] <- yy[fit$data[[by]] == by_levels[b]] -
            predict(fit, newdata = fit$data[fit$data[[by]] == by_levels[b],]) +
            approx(new_data[[coef_name]], preds[[b]],
                   xout = fit$data[[coef_name]][fit$data[[by]] == by_levels[b]])$y
        }


        # Validate plot.args
        if(!length(plot.args$ylab))
          plot.args$ylab <- outcome_name

        if(!length(plot.args$xlab))
          plot.args$xlab <- coef_name

        if(!length(plot.args$ylim))
          plot.args$ylim <- range(preds, resid, na.rm = TRUE)

        if(!length(plot.args$xlim))
          plot.args$xlim <- range(new_data[[coef_name]])

        cols <- pal(length(by_levels))
        plot(new_data[[coef_name]], preds[[1]], type = "n", lwd = 2,
             xlab = plot.args$xlab, ylab = plot.args$ylab,
             ylim = plot.args$ylim, xlim = plot.args$xlim, ...)



        lty <- rep(1, length(by_levels))
        for(i in 1:length(by_levels)) {
          ## Check for overlapping lines
          if(i > 2) {
            overlap <- which(sapply(1:(i-1), function(j) isTRUE(all.equal(preds[[j]], preds[[i]]))))
          } else
            overlap <- NULL

          if(length(overlap))
            lty[i] = 2

          lines(new_data[[coef_name]], preds[[i]], lwd = 2, col = cols[i], lty = lty[i])
        }

        points(fit$data[[coef_name]], unlist(resid), pch = 20, col = pt_cols)


        leg_lab <- paste0(by, ": ", unname(by_levels))
        legend(legend_location, legend = leg_lab, col = cols, lwd = 2, bty = "n", lty = lty)

        ## Factor main covariate
      } else
        stop("Cannot yet handle factor effect for both main covariate and 'by' covariate")

      ## Numeric by variable
    } else if(inherits(fit$data[[by]], what = "numeric")) {


      ## IF it's a numeric variable
      if(!cfact) {
        if(missing(by_levels)) {
          by_levels <- quantile(fit$data[[by]], c(.1, .5, .9), na.rm = TRUE)
          closest_level <- apply(sapply(by_levels, function(x) (abs(x - fit$data[[by]]))), 1, which.min)
        }

        resid <- preds <- list()
        for(b in 1:length(by_levels)) {
          new_data[[by]] <- rep(by_levels[[b]], nn)

          preds[[b]] <- predict(fit, newdata=new_data)
          resid[[b]] <- yy[closest_level == b] -
            predict(fit, newdata = fit$data[closest_level == b,]) +
            approx(new_data[[coef_name]], preds[[b]], xout = fit$data[[coef_name]][closest_level == b])$y
        }


        # Validate plot.args
        if(!length(plot.args$ylab))
          plot.args$ylab <- outcome_name

        if(!length(plot.args$xlab))
          plot.args$xlab <- coef_name

        if(!length(plot.args$ylim))
          plot.args$ylim <- range(preds, na.rm = TRUE)

        if(!length(plot.args$xlim))
          plot.args$xlim <- range(new_data[[coef_name]])

        cols <- pal(length(by_levels))
        plot(new_data[[coef_name]], preds[[1]], type = "n", lwd = 2,
             xlab = plot.args$xlab, ylab = plot.args$ylab,
             ylim = plot.args$ylim, xlim = plot.args$xlim, ...)

        ll <- lapply(1:length(by_levels), function(i) {
          lines(new_data[[coef_name]], preds[[i]], lwd = 2, col = cols[i])
        })

        pp <- lapply(1:length(by_levels), function(i) {
          points(fit$data[[coef_name]][closest_level == i], resid[[i]], pch = 20, col = cols[i])
        })

        leg_lab <- paste0(by, ": ", round(unname(by_levels),2))
        legend(legend_location, legend = leg_lab, col = cols, lwd = 2, bty = "n")

        ## factor covariate by numeric one
      } else {
        stop("Cannot yet handle factor covariate by numeric covariate")
      }
    }
  }
  return(invisible(NULL))
}
