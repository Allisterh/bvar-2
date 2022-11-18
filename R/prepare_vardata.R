# This is an auxiliary function.
# It sets the data in SUR form
# All data should be submitted in data frame form
# The function inputs include the data in a data frame form, the number of lags and logical values about the inclusion of constant, trend or quadratic trend.
# The inclusion of Exogenous variables is also considered
# The function returns the RHS and teh LHS of a VAR equation and some of the arguments (final number of observations,
# the number of the exogenous variables, the number of the deterministic terms
VAR_data_prep = function(data , lags = 1, const = T, trend = F, trend_qua = F, ex = NULL, ex_lag = NULL, xlags = 0 ){

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
  ndet = 0 # Sum the number of the deterministic terms
  nexo = 0 # Sum the number of the exogenous variables
  n = ncol(data)  # number of variables
  t = nrow(data)  # number of observations
  dt = NULL

  # Prepare the data
  # separate exogenous from endogenous variables and calculate the lags
  y = data

  # calculate lags for endogenous
  x = dplyr::lag(y, 1)
  names(x) = str_c(names(x),1, sep = "_")
  if (lags >1 ){
    for (i in 2:lags) {
      lag_t = dplyr::lag(y, i)
      names(lag_t) = str_c(names(lag_t),i, sep = "_")
      x = cbind(x,  lag_t)
    }
  }

  # check for constant
  if (const == T) {
    dt = cbind(c = rep(1, t)) # First column filled with 1s
    ndet = ndet + 1
  }

  # check for trend
  if (trend == T) {
    dt = cbind(dt, "trend" = 1:nrow(y)) # First column filled with 1s
    ndet = ndet + 1
  }

  # check for quadratic trend
  if (trend_qua == T) {
    dt = cbind( dt, "q_trend"  = c(1:nrow(y))^2) # First column filled with 1s
    ndet = ndet + 1
  }

  # calculate lags for the exogenous
  if(!is.null(ex_lag)){

    te_ex_lag = ex_lag
    if (xlags > 0){
      for (i in 1:xlags) {
        xlag_t = dplyr::lag(te_ex_lag, i)
        names(xlag_t) = str_c(names(xlag_t),i, sep = "_")
        ex_lag = cbind(ex_lag,  xlag_t)
      }
    }
    nexo = nexo + ncol(ex_lag)

    # bind the exogenous variables with the lagged matrix of endogenous
    x = cbind(ex_lag,x)
  }

  if(!is.null(ex)){
    x = cbind(ex,x)
    nexo = nexo + ncol(ex)
  }

  if(!is.null(dt)){
    x = cbind(dt,x)
  }

  y_t =  as.matrix(y[(lags+1):t,]) # LHS
  x_t =  as.matrix(x[(lags+1):t,]) # RHS

  nexo = ndet + nexo
  names_endo  = names(y)
  L_R_HS = list(y_t,x_t,n,nexo = nexo, ndet, t, t_lag = nrow(y_t), names_endo)
  names(L_R_HS) = c("LHS", "RHS","n", "nexo","ndet","t_obs", "t_lag", "names_endo")
  return(L_R_HS)
}
