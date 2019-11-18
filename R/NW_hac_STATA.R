#' NW_hac_STATA
#' tbd
#'
#' @param vars tbd
#' @param lags tbd
#'
#' @return Sigma: tbd
#'
#' @export
#'

NW_hac_STATA <- function(vars, lags) {
  Sigma0 = (1 / (dim(vars)[1])) * (t(vars) %*% vars)
  Sigma = Sigma0
  if (lags >= 1) {
    for (n in 1:lags) {
      Sigma = Sigma + (1 - n / (lags + 1)) * (Sigma_cov(vars, n) + t(Sigma_cov(vars, n)))
    }
  }
  return(Sigma)
}

Sigma_cov <- function(vars, k) {
  return((1 / (dim(vars)[1])) * t(vars[1:(nrow(vars) - k), ]) %*% vars[(1 +
                                                                          k):nrow(vars), ])
}


