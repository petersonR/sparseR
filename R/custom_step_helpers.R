## Custom Stepwise helper functions (
# these are helper functions for the custom step functions, which allows MASS::stepAIC to work with other penalties)

## Some internal functions needed
safe_pchisq <- function (q, df, ...){
  df[df <= 0] <- NA
  pchisq(q = q, df = df, ...)
}

safe_pf <- function (q, df1, ...) {
  df1[df1 <= 0] <- NA
  pf(q = q, df1 = df1, ...)
}

run_steps <- function (object, scope, scale = 0, direction = c("both", "backward", "forward"),
                       trace = 1, keep = NULL, steps = 1000, use.start = FALSE, k = 2, penalty = NULL, ...) {

  if(!length(penalty)) {
    penaltyfn <- extractAIC
  } else if(cid <- "character" %in% class(penalty)) {
    p_char <- penalty[[cid]]
    penaltyfn <- ifelse(p_char == "AIC", extractAIC,
                        "Default not recognized")
  } else if("function" %in% class(penalty)) {
    penaltyfn <- penalty
  } else stop("penalty not recognized")

  mydeviance <- function(x, ...) {
    dev <- deviance(x)
    if (!is.null(dev))
      dev
    else
      extractAIC(x, k = 0)[2L]
  }

  cut.string <- function(string) {
    if (length(string) > 1L)
      string[-1L] <- paste("\n", string[-1L], sep = "")
    string
  }
  re.arrange <- function(keep) {
    namr <- names(k1 <- keep[[1L]])
    namc <- names(keep)
    nc <- length(keep)
    nr <- length(k1)
    array(unlist(keep, recursive = FALSE), c(nr, nc), list(namr,
                                                           namc))
  }
  step.results <- function(models, fit, object, usingCp = FALSE) {
    change <- sapply(models, "[[", "change")
    rd <- sapply(models, "[[", "deviance")
    dd <- c(NA, abs(diff(rd)))
    rdf <- sapply(models, "[[", "df.resid")
    ddf <- c(NA, abs(diff(rdf)))
    AIC <- sapply(models, "[[", "AIC")
    heading <-
      c(
        "Stepwise Model Path \nAnalysis of Deviance Table",
        "\nInitial Model:",
        deparse(formula(object)),
        "\nFinal Model:",
        deparse(formula(fit)),
        "\n"
      )
    aod <-
      data.frame(
        Step = change,
        Df = ddf,
        Deviance = dd,
        `Resid. Df` = rdf,
        `Resid. Dev` = rd,
        IC = IC,
        check.names = FALSE
      )
    attr(aod, "heading") <- heading
    class(aod) <- c("Anova", "data.frame")
    fit$anova <- aod
    fit
  }
  Terms <- terms(object)
  object$formula <- Terms
  if (inherits(object, "lme"))
    object$call$fixed <- Terms
  else if (inherits(object, "gls"))
    object$call$model <- Terms
  else
    object$call$formula <- Terms
  if (use.start)
    warning("'use.start' cannot be used with R's version of 'glm'")
  md <- missing(direction)
  if(direction != "forward")
    stop("run_steps can currently only use direction = 'forward'")

  forward <- TRUE

  if (missing(scope)) {
    fadd <- attr(Terms, "factors")
  } else {
    if (is.list(scope)) {
      fadd <- if (!is.null(fadd <- scope$upper))
        attr(terms(update.formula(object, fadd)), "factors")
    } else {
      fadd <- if (!is.null(fadd <- scope))
        attr(terms(update.formula(object, scope)), "factors")
    }
  }
  models <- vector("list", steps)
  if (!is.null(keep))
    keep.list <- vector("list", steps)
  n <- nobs(object, use.fallback = TRUE)
  fit <- object
  bIC <- penaltyfn(fit, ...)

  edf <- bIC[1L]
  bIC <- bIC[2L]
  if (is.na(bIC))
    stop("IC is not defined for this model, cannot proceed")
  if (bIC == -Inf)
    stop("IC is -infinity for this model, cannot proceed")
  nm <- 1
  Terms <- terms(fit)
  if (trace) {
    cat("Start:  IC=",
        format(round(bIC, 2)),
        "\n",
        cut.string(deparse(formula(fit))),
        "\n\n",
        sep = "")
    utils::flush.console()
  }
  models[[nm]] <- list(
    deviance = mydeviance(fit),
    df.resid = n -
      edf,
    change = "",
    IC = bIC
  )
  if (!is.null(keep))
    keep.list[[nm]] <- keep(fit, bIC)
  usingCp <- FALSE
  while (steps > 0) {
    steps <- steps - 1
    IC <- bIC
    ffac <- attr(Terms, "factors")
    if (!is.null(sp <-
                 attr(Terms, "specials")) && !is.null(st <- sp$strata))
      ffac <- ffac[-st,]
    scope <- factor.scope(ffac, list(add = fadd))
    aod <- NULL
    change <- NULL
    if (is.null(change)) {
      if (forward && length(scope$add)) {
        aodf <- my_addterm(
          object = fit,
          scope$add,
          scale = scale,
          trace = max(0, trace - 1),
          k = k,
          penaltyfn = penaltyfn,
          ...
        )
        rn <- row.names(aodf)
        row.names(aodf) <- c(rn[1L], paste("+", rn[-1L],
                                           sep = " "))
        if (is.null(aod)) {
          aod <- aodf
        } else aod <- rbind(aod, aodf[-1, , drop = FALSE])
      }
      attr(aod, "heading") <- NULL
      if (is.null(aod) || ncol(aod) == 0)
        break
      nzdf <- if (!is.null(aod$Df))
        aod$Df != 0 | is.na(aod$Df)
      aod <- aod[nzdf,]
      if (is.null(aod) || ncol(aod) == 0)
        break
      nc <- match(c("Cp", "IC"), names(aod))
      nc <- nc[!is.na(nc)][1L]
      o <- order(aod[, nc])
      if (trace) {
        print(aod[o,])
        utils::flush.console()
      }
      if (o[1L] == 1)
        break
      change <- rownames(aod)[o[1L]]
    }
    usingCp <- match("Cp", names(aod), 0) > 0
    fit <- update(fit, paste("~ .", change), evaluate = FALSE)
    fit <- eval.parent(fit)
    nnew <- nobs(fit, use.fallback = TRUE)
    if (all(is.finite(c(n, nnew))) && nnew != n)
      stop("number of rows in use has changed: remove missing values?")
    Terms <- terms(fit)
    bIC <- penaltyfn(fit, ...)
    edf <- bIC[1L]
    bIC <- bIC[2L]
    if (trace) {
      cat("\nStep:  IC=",
          format(round(bIC, 2)),
          "\n",
          cut.string(deparse(formula(fit))),
          "\n\n",
          sep = "")
      utils::flush.console()
    }
    if (bIC >= IC + 1e-07)
      break
    nm <- nm + 1
    models[[nm]] <- list(
      deviance = mydeviance(fit),
      df.resid = n -edf,
      change = change,
      IC = bIC
    )
    if (!is.null(keep))
      keep.list[[nm]] <- keep(fit, bIC)
  }
  if (!is.null(keep))
    fit$keep <- re.arrange(keep.list[seq(nm)])
  step.results(models = models[seq(nm)], fit, object, usingCp)
}


fitall <- function(y, X, method = "lm", ...) {
  data <- cbind(y=y, X)
  X <- as.data.frame(X)

  combs <- do.call(expand.grid, rep(list(c(FALSE, TRUE)), ncol(X)))

  vars <- apply(combs, 1, function(i) names(X)[i])
  vars[[1]] <- "1"
  form <- paste("y ~ ", lapply(vars, paste, collapse=" + "), sep = "")
  form <- lapply(form, as.formula)

  method <- as.name(method)
  fitmodel <- function(f) {
    eval(substitute(method(f, data = data, model = FALSE, ...),
                    list(f = f, method = method)))
  }

  models <- lapply(form, fitmodel)
  models
}

# See stats:::check_exact()
check_exact <- function (object) {
  w <- object$weights
  if (is.null(w)) {
    mss <- sum(object$fitted.values^2)
    rss <- sum(object$residuals^2)
  }
  else {
    mss <- sum(w * object$fitted.values^2)
    rss <- sum(w * object$residuals^2)
  }
  if (rss < 1e-10 * mss)
    warning("attempting model selection on an essentially perfect fit is nonsense",
            call. = FALSE)
}

# See stats:::deviance.lm()
deviance.lm <- function (object, ...)
  sum(weighted.residuals(object)^2, na.rm = TRUE)


my_add1 <- function (object, scope, scale = 0,
                     test = c("none"),
                     x = NULL, k = 2, penaltyfn, ...) {

  check_exact(object)
  if (missing(scope) || is.null(scope))
    stop("no terms in scope")
  if (!is.character(scope))
    scope <- add.scope(object, update.formula(object, scope))
  if (!length(scope))
    stop("no terms in scope for adding to object")
  oTerms <- attr(object$terms, "term.labels")
  int <- attr(object$terms, "intercept")
  ns <- length(scope)
  y <- object$residuals + object$fitted.values
  dfs <- numeric(ns + 1)
  RSS <- numeric(ns + 1)
  ic <- numeric(ns + 1)
  names(dfs) <- names(RSS) <- names(ic) <- c("<none>", scope)
  add.rhs <- paste(scope, collapse = "+")
  add.rhs <- eval(parse(text = paste("~ . +", add.rhs), keep.source = FALSE))
  new.form <- update.formula(object, add.rhs)
  Terms <- terms(new.form)
  if (is.null(x)) {
    fc <- object$call
    fc$formula <- Terms
    fob <- list(call = fc, terms = Terms)
    class(fob) <- oldClass(object)
    m <- model.frame(fob, xlev = object$xlevels)
    x <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
    offset <- model.offset(m)
    wt <- model.weights(m)
    oldn <- length(y)
    y <- model.response(m, "numeric")
    newn <- length(y)
    if (newn < oldn)
      warning(sprintf(ngettext(newn, "using the %d/%d row from a combined fit",
                               "using the %d/%d rows from a combined fit"),
                      newn, oldn), domain = NA)
  }
  else {
    wt <- object$weights
    offset <- object$offset
  }
  n <- nrow(x)
  Terms <- attr(Terms, "term.labels")
  asgn <- attr(x, "assign")
  ousex <- match(asgn, match(oTerms, Terms), 0L) > 0L
  if (int)
    ousex[1L] <- TRUE
  iswt <- !is.null(wt)
  X <- x[, ousex, drop = FALSE]
  if (iswt) {
    z <- lm.wfit(X, y, wt, offset = offset)
  }  else z <- lm.fit(X, y, offset = offset)
  dfs[1L] <- z$rank
  z$nobs <- length(y)
  class(z) <- "lm"
  RSS[1L] <- deviance(z)
  ic[1L] <- penaltyfn(z, ...)[2]

  sTerms <- sapply(strsplit(Terms, ":", fixed = TRUE), function(x) paste(sort(x),
                                                                         collapse = ":"))
  for (tt in scope) {
    stt <- paste(sort(strsplit(tt, ":")[[1L]]), collapse = ":")
    usex <- match(asgn, match(stt, sTerms), 0L) > 0L
    X <- x[, usex | ousex, drop = FALSE]
    z <- if (iswt)
      lm.wfit(X, y, wt, offset = offset)
    else lm.fit(X, y, offset = offset)
    class(z) <- "lm"
    z$nobs <- length(y)
    dfs[tt] <- z$rank
    RSS[tt] <- deviance(z)
    ic[tt] <- penaltyfn(z, ...)[2]
  }

  dfs <- dfs - dfs[1L]
  dfs[1L] <- NA
  aod <- data.frame(Df = dfs, `Sum of Sq` = c(NA, RSS[1L] - RSS[-1L]),
                    RSS = RSS, IC = ic, row.names = names(dfs),
                    check.names = FALSE)
  if (scale > 0)
    names(aod) <- c("Df", "Sum of Sq", "RSS", "Cp")
  test <- match.arg(test)
  head <- c("Single term additions", "\nModel:", deparse(formula(object)),
            if (scale > 0) paste("\nscale: ", format(scale), "\n"))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}

my_addterm <- function (object, scope, scale = 0,
                        test = c("none"),
                        k = 2, sorted = FALSE,
                        penaltyfn = penaltyfn,
                        ...) {

  if (missing(scope) || is.null(scope))
    stop("no terms in scope")
  aod <- my_add1(object, scope = scope, scale = scale, k=k, penaltyfn = penaltyfn, ...)
  dfs <- c(0, aod$Df[-1L]) + object$rank
  RSS <- aod$RSS
  n <- length(object$residuals)

  o <- if (sorted)
    order(aod$IC)
  else seq_along(aod$IC)

  if (scale > 0)
    names(aod) <- c("Df", "Sum of Sq", "RSS", "Cp")
  test <- match.arg(test)
  aod <- aod[o, ]
  head <- c("Single term additions", "\nModel:", deparse(formula(object)))
  if (scale > 0)
    head <- c(head, paste("\nscale: ", format(scale), "\n"))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}
