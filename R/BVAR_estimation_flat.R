


BVAR_estimation_flat = function(data, lags = 4, draws, const = T, trend = F, trend_qua = F, ex = NULL, ex_lag = NULL, xlags = 0) {
  data = df[,1:2]; lags =4; draws = 10;const= T;trend =T; trend_qua = T;  ex = NULL ; ex_lag = df[,3, drop = F]; xlags = 4; i=1
  vardata = var_data_prepare(data, lags = lags, const = const, trend = trend, trend_qua = trend_qua , ex = ex, ex_lag = ex_lag, xlags = xlags)
  # Take the arguments from the
  n = vardata$number_of_variables
  nexo = vardata$number_of_exogenous
  t_lag = vardata$number_of_obs_lags

  CM = array(0, list(n * lags , n * lags, draws))
  Sigma = array(0, list(n, n, draws))
  res = array(0, list(t_lag, n, draws))
  fit = array(0, list(t_lag, n, draws))
  beta_s = array(0, list((n*lags + nexo),n, draws), dimnames=list( colnames(vardata$RHS) , vardata$names_endo , NULL ))
  # Reduced-form VAR estimation:

  for (i in 1:draws) {
    stability  = F

    while (!stability) {
      # BRAR_draw draws from normal & inverse wishart distribution
      bvar_hat =  BVAR_draw(data = data,lags = lags, const = const, trend = trend, trend_qua = trend_qua, ex = ex, ex_lag = ex_lag, xlags = xlags)
      # save the returned values of the function
      beta_s[, , i] = bvar_hat$beta
      CM[, , i] = bvar_hat$CM
      Sigma[, , i] = bvar_hat$Sigma
      res[, , i] =  bvar_hat$res
      fit[, , i] = bvar_hat$y_fit
      stability = bvar_hat$stability
    }

  }

  bvar = list(beta_s, Sigma, CM, res, fit, vardata, draws)
  names(bvar) = c("posterior","Sigma", "CM", "res", "yfit", "vardata", "draws")
  return(bvar)

}

