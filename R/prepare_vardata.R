# This is an auxiliary function.
# It sets the data in SUR form
# All data should be submitted in data frame form
# The function inputs include the data in a data frame form, the number of lags and logical values about the inclusion of constant, trend or quadratic trend.
# The inclusion of Exogenous variables is also considered
# The function returns the RHS and teh LHS of a VAR equation and some of the arguments (final number of observations,
# the number of the exogenous variables, the number of the deterministic terms
#' Title
#'
#' @param data
#' @param lags
#' @param const
#' @param trend
#' @param trend_qua
#' @param ex
#' @param ex_lag
#' @param xlags
#'
#' @return
#' @export
#'
#' @examples
#'

var_data_prepare = function(data , lags = 1, const = T, trend = F, trend_qua = F, ex = NULL, ex_lag = NULL, xlags = 0 ){

  # Check the consistency of the inputs
  # Test the data input
  if (!is.data.frame(data) & !is.tbl(data)) {
    stop("Please provide the argument 'data' either as a 'data.frame' or a 'tbl' object.")
  }

  if (any(is.na(data))) {
    stop("The argument data contains NAs. Function cannot handle mssing values.")
  }

  if (!is.numeric(as.matrix(data))) {
    stop("The argument data contains non numeric values. Please check you data.")
  }

  # Test the lags input
  if ( !is.numeric(lags)) {
    stop("Object lags must be a positive integer greater or equal to 1")
  }

  if (!is_scalar_atomic(lags) | !lags%%1 == 0 | lags < 1) {
    stop("Object lags must be a positive integer greater or equal to 1")
  }

  # Test the  inclusion of deterministic terms
  if( !is.logical(const) ) {
    stop("argument const takes a logical value. The default value is 'TRUE'")
  }

  if( !is.logical(trend)) {
    stop("argument trend takes a logical value. The default value is 'FALSE'")
  }

  if( !is.logical(trend_qua)) {
    stop("trend_qua takes a logical value. The default value is 'FALSE'")
  }

  # Test the ex input
  if (!is.data.frame(ex) & !is.tbl(ex) & !is.null(ex)) {
    stop("Please provide the argument 'ex' either as a 'data.frame' or a 'tbl' object or a matrix order. The default value is 'NULL'")
  }

  if (!is.null(ex)) {
    if (!is.numeric(as.matrix(ex))) {
      stop("The argument 'ex' does not include numeric values")
    }

    if (any(is.na(ex))) {
      stop("The argument 'ex' contains NAs. Function cannot handle mssing values.")
    }
  }

  # Test the ex_lag input
  if (!is.data.frame(ex_lag) & !is.tbl(ex_lag) & !is.null(ex_lag)) {
    stop("Please provide the argument 'ex_lag' either as a 'data.frame' or a 'tbl' object or a matrix order.")
  }

  if (!is.null(ex_lag)) {
    if (!is.numeric(as.matrix(ex_lag))) {
      stop("The argument 'ex_lag' does not include numeric values")
    }

    if (any(is.na(ex_lag))) {
      stop("The argument 'ex_lag' contains NAs. Function cannot handle mssing values.")
    }
  }

  # Test the compatibility of the objects
  if (!is.null(ex)){
    if ( !dim(data)[1] == dim(ex)[1] ){
      stop("Arguments 'data' and 'ex' are not of the same langht")
    }
  }

  if (!is.null(ex_lag)){
    if ( !dim(data)[1] == dim(ex_lag)[1] ){
      stop("Arument 'data' and 'ex_lag' are not of the same langht")
    }
  }

  # Test the xlags input
  if ( !is.numeric(xlags) | !is_scalar_atomic(xlags) | !xlags%%1 == 0 | xlags < 0 ) {
    stop("Object xlags must be a positive integer greater or equal to 1")
  }
  if (xlags > lags) {
    stop("Argument 'xlags' must be lower than argument 'lag'.")
  }

  # Start data preparation #
  number_of_deterministic <- 0 # Sum the number of the deterministic terms

  number_of_exogenous     <- 0 # Sum the number of the exogenous variables

  number_of_variables     <- ncol(data)  # number of variables

  number_of_observations  <- nrow(data)  # number of observations

  # Initialize the block of the deterministc terms
  deterministic_block     <- NULL

  # separate exogenous from endogenous variables and calculate the lags
  y = data

  # calculate lags for endogenous
  x  <-  dplyr::lag(y, 1)

  names(x) <- paste(names(x), 1, sep = "_")

  if ( lags >1 ) {

    for (i in 2:lags) {

      lag_temp <- dplyr::lag(y, i)

      names(lag_temp) <- paste(names(lag_temp),i, sep = "_")

      x <- cbind(x,  lag_temp)

    }

  }

  # check for constant
  if ( const == T ) {

    deterministic_block <- cbind(c = rep(1, number_of_observations)) # First column filled with 1s

    number_of_deterministic <- number_of_deterministic + 1

  }

  # check for trend
  if ( trend == T ) {

    deterministic_block <- cbind(deterministic_block, "trend" = 1:nrow(y)) # Second column filled with linear trend

    number_of_deterministic <- number_of_deterministic + 1

  }

  # check for quadratic trend
  if ( trend_qua == T ) {

    deterministic_block = cbind(deterministic_block, "q_trend"  = c(1:nrow(y))^2) # third column filled with quadratic trend

    number_of_deterministic <- number_of_deterministic + 1

  }

  # calculate lags for the exogenous
  if( !is.null(ex_lag) ) {

    temp_ex_lag = ex_lag

    if ( xlags > 0 ) {

      for (i in 1:xlags) {

        xlag_temp = dplyr::lag(temp_ex_lag, i)

        names(xlag_temp) = str_c(names(xlag_temp),i, sep = "_")

        ex_lag = cbind(ex_lag,  xlag_temp)

      }
    }
    number_of_exogenous = number_of_exogenous + ncol(ex_lag)

    # bind the exogenous variables with the lagged matrix of endogenous
    x = cbind(ex_lag,x)

  }

  if( !is.null(ex) ) {

    x = cbind(ex,x)

    number_of_exogenous = number_of_exogenous + ncol(ex)

  }


  if( !is.null(deterministic_block) ){

    x = cbind(deterministic_block,x)

  }

  # I turn the data into matrices to facilitate the step of the calculations
  y_lhs =  as.matrix(y[(lags+1):number_of_observations,]) # LHS
  x_rhs =  as.matrix(x[(lags+1):number_of_observations,]) # RHS

  # Total number of exogenous variables
  number_of_exogenous = number_of_deterministic + number_of_exogenous

  #
  names_endogenous  = names(y)

  # Function will return a list with the VAR data and a number of parameters
  # that will be used in the calculations
  vardata = list(
    "y_lhs" = y_lhs,
    "x_lhs" = x_rhs,
    "number_of_variables"     = number_of_variables,
    "number_of_exogenous"     = number_of_exogenous,
    "number_of_deterministic" = number_of_deterministic,
    "number_of_observations"  = number_of_observations,
    "number_of_obs_lags"      = nrow(y_lhs),
    "names_of_endogenous"     = names_endogenous
  )


  return(vardata)
}
