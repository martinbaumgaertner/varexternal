#' MARep
#'
#' Transforms the A(L) parameters of a reduced-form VAR into the coefficients C of the MA representation.
#'
#' @param AL VAR model coefficients
#' @param p lag order
#' @param hori forecast horizon
#'
#' @return C MA representation coefficients
#'
#' @export
MARep <- function(AL, p, hori) {
  n         = dim(AL)[1]

  vecAL     = array(AL, c(n, n, p))

  vecALrevT = array(rep(matrix(0,n,n), hori), c(n, n, hori))
  for (ihori in 1:hori) {
    if (ihori < (hori - p) + 1) {
      vecALrevT[, , ihori] = matrix(0,n,n)
    } else{
      vecALrevT[, , ihori] = t(vecAL[, , (hori - ihori) + 1])
    }

  }
  vecALrevT     = array(vecALrevT, c(n, n * hori))

  C             = repmat(vecAL[, , 1], 1, hori)

  for (ihori in 1:(hori - 1)) {
    C[, ((n * ihori) + 1):(n * (ihori + 1))] = cbind(diag(n), C[, 1:(n * ihori)]) %*% t(vecALrevT[, ((hori *
                                                                                                       n - (n * (ihori + 1))) + 1):ncol(vecALrevT)])

  }
  return(C)
}


