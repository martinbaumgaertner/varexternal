#' CovAhat_Sigmahat_Gamma
#'
#' Estimates the asymptotic covariance matrix
#'
#' @param p Number of lags in the VAR model
#' @param X VAR "right-hand" variables
#' @param Z external instrument
#' @param eta eta
#' @param lags Newey-West lags
#'
#' @return WHataux asymptotic variance of
#' @return WHat asymptotic variance of
#' @return V Matrix such that
#'
#'
#' @export

CovAhat_Sigmahat_Gamma<-function(p,X,Z,eta,lags){
  n      = nrow(eta)
  k      = 1 # ncol(Z)
  m      = ncol(X) - (n * p)
  XSVARp = X
  matagg = t(cbind(XSVARp, t(eta), Z))
  T1aux  = ncol(eta)
  T2aux  = nrow(matagg)
  etaaux = array(eta, c(n, 1, T1aux))
  mataggaux = aperm(array(matagg, c(T2aux, 1, T1aux)), c(2, 1, 3))

  d <-array(data = NA,
          dim = c(dim(etaaux)[1], dim(mataggaux)[2], dim(mataggaux)[3]),
          dimnames = NULL)
  daux<-outer(etaaux, mataggaux, FUN = "*")
  for (i in 1:dim(mataggaux)[3]) {
    d[, , i] <- daux[, , i, , , i]
  }
  d1 <- array(data = NA, dim = c(dim(etaaux)[1], dim(mataggaux)[2]))
  for (i in 1:dim(etaaux)[1]) {
    for (j in 1:dim(mataggaux)[2]) {
      d1[i, j] <- mean(d[i, j, ])
    }
  }
  auxeta <- sweep(d, 1:2, d1, FUN = "-")

  vecAss1 = array(auxeta, c((n * m) + (p * (n ^ 2)) + n ^ 2 + (n * k), 1, T1aux))

  AuxHAC1  = vecAss1[1:dim(vecAss1)[1], , ]

  AuxHAC2  = t(array(AuxHAC1, c(dim(vecAss1)[1], dim(vecAss1)[3])))

  AuxHAC3  = NW_hac_STATA(AuxHAC2, lags)

  WhatAss1 = AuxHAC3

  I = diag(n)
  V = t(kronecker(I[1, ], I))

  for (i_vars in 2:n) {
    d = i_vars - 1
    a = I[i_vars, ]
    V    = rbind(V, t(kronecker(a, I))[-c(1:d), ])
  }


  Q1       = (t(XSVARp) %*% XSVARp / T1aux)

  Q2       = t(Z) %*% XSVARp / T1aux


  Shat = rbind(cbind(kronecker(cbind(
    matrix(0,n * p,m), diag(n * p)
  ) %*% (solve(
    Q1
  )), diag(n)), matrix(0,(n ^ 2) * p, n ^ 2 + (k * n))),
  cbind(matrix(0,n * (n + 1) / 2, ((
    n ^ 2
  ) * p) + (n * m)), V, matrix(0,n * (n + 1) / 2, k * n)),
  cbind(-kronecker(Q2 %*% (solve(
    Q1
  )), diag(n)), matrix(0,k * n, n ^ 2), diag(k * n)))

  WHataux  = (Shat) %*% (WhatAss1) %*% (t(Shat))

  WHat <- rbind(cbind(WHataux[1:((n ^ 2) * p), 1:((n ^ 2) * p)],
                      WHataux[1:((n ^ 2) * p), (((n ^ 2) * p) + (n * (n +
                                                                        1) / 2) + 1):ncol(WHataux)]),
                cbind(t(WHataux[1:((n ^ 2) * p), (((n ^ 2) * p) + (n * (n +
                                                                          1) / 2) + 1):ncol(WHataux)]),
                      WHataux[(((n ^ 2) * p) + (n * (n + 1) / 2) + 1):nrow(WHataux), (((n ^
                                                                                          2) * p) + (n * (n + 1) / 2) + 1):ncol(WHataux)]))

  return(list(WHataux = WHataux,
              WHat = WHat,
              V = V))
}





