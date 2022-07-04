#' Helper function to help set up penalties
#' @param varnames names of the covariates in the model matrix
#' @param poly max polynomial considered
#' @param poly_prefix what comes before the polynomial specification in these
#'   varnames?
#' @param int_sep What denotes the multiplication for interactions?
#' @param pool Should polynomials and interactions be pooled?
#' @param gamma How much should the penalty increase with group size (0.5
#'   assumes equal contribution of prior information)
#' @param cumulative_k Should penalties be increased cumulatively as order
#'   interaction increases? (only used if !pool)
#' @param cumulative_poly Should penalties be increased cumulatively as order
#'   polynomial increases? (only used if !pool)
#' @details This is primarily a helper function for sparseR, but it may be
#'   useful if doing the model matrix set up by hand.
#'
#' @return a list of relevant information for the variables, including:
#'
#' \item{penalties}{the numeric value of the penalties}
#' \item{vartype}{Variable type (main effect, order k interaction, etc)}
#' \item{varname}{names of variables}
#'
#' @export
get_penalties <- function(varnames, poly, poly_prefix = "poly_", int_sep = "\\:",
                          pool = FALSE, gamma = .5, cumulative_k = FALSE,
                          cumulative_poly = TRUE) {

  ## Get count of interaction seps for each varname
  kk <- sapply(gregexpr(int_sep, varnames), function(x) {
    if(x[1] == -1)
      return(NA)
    length(x)
  })

  ## Get polynomial group for each varname
  p_id <- gregexpr(poly_prefix, varnames)
  pp <- sapply(1:length(p_id), function(i) {
    x <- p_id[[i]]
    if(x[1] == -1)
      return(NA)

    ## Requires accomodation for when polynomials are greater than 1 digit
    nchar <- nchar(gsub("\\\\", "", poly_prefix))
    as.numeric(
      substr(varnames[i], x + nchar,
             x + nchar + (poly > 9) + (poly > 99)))
  })

  # Code dummy effects as 0, if necessary
  if(any(idx <- is.na(kk) & is.na(pp))) {
    pp[idx] <- 1
    kk[idx] <- 0
  }

  ## Code polynomial of order 1 as k = 0
  if(any(idx <- pp == 1)) {
    kk[idx] <- 0
  }


  if(!pool) {
    ## Calculate totals and respective penalties for interactions
    if(cumulative_k) {
      kpen <- cumsum(table(kk)) ^ gamma
    } else
      kpen <- (table(kk)) ^ gamma

    # Expand k penalties for all variables
    kkpen <- as.numeric(as.character(factor(kk, levels = names(kpen), labels = kpen)))

    ## Set NA values (polynomial of order > 1) to the same penalty as main effects (for now)
    kkpen[is.na(kkpen)] <- kkpen[which(kk==0)[1]]

    ## Get factor increase for poly vars
    if(cumulative_poly) {
      ppen <- (cumsum(table(pp))) ^ gamma
      ppen <- ppen / min(ppen)

    } else {
      ppen <- (table(pp)) ^ gamma
      ppen <- ppen / min(ppen)
    }

    # Expand p penalties for all variables
    pppen <- as.numeric(as.character(factor(pp, levels = names(ppen), labels = ppen)))
    pppen[is.na(pppen)] <- 1

    pen <- kkpen * pppen

  } else {

    ## Create pooled values
    pool_pen <- pp
    pool_pen[is.na(pp)] <- kk[is.na(pp)]+1

    # Get factor increase for poly vars
    if(cumulative_poly) {
      ppen <- (cumsum(table(pool_pen))) ^ gamma
    } else {
      ppen <- (table(pool_pen)) ^ gamma
    }

    # Expand penalties for all variables
    pen <- as.numeric(as.character(factor(pool_pen, levels = names(ppen), labels = ppen)))
  }

  vartypes <- rep("Main effect", length(pen))
  vartypes <- ifelse(kk > 0, paste("Order", kk, "interaction"), vartypes)
  vartypes <- ifelse(!is.na(pp) & pp > 1, paste("Order", pp, "polynomial"), vartypes)

  list(penalties = pen, vartypes = vartypes, varnames = varnames)
}
