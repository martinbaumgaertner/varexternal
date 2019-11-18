#' RForm_VAR
#'
#' Provides reduced form estimators of a VAR(p) model
#'
#' @param TSL matrix of time series
#' @param p number of lags in the VAR model
#' @param W Matrix of exogenous regressors
#'
#' @return AL: Least-squares estimator of the VAR coefficients
#' @return Sigma: Least-squares estimator of the VAR residuals
#' @return eta: VAR model residuals
#' @return X: Matrix of VAR covariates
#' @return Y: VAR matrix of endogenous regressors
#'
#' If no exogenous regressors are specified, our estimation always includes a constant.
#'
#' @seealso https://github.com/jm4474/SVARIV
#' @seealso https://uc4384f22718973d3d912d614f2e.dl.dropboxusercontent.com/cd/0/inline2/Aslwi_FA9z6v37dfxMeUP3mEaP-a4ZA8RFrXd5h_EymM1lqGy-UUwGJB3TxqRxluMhE9ebnvHEMZktfcUmJ1v4mGSxeS_ZEDJBCrjGvyopAnP44ClPLJ6kmM36F557kHnqKQjn___0XTNvPaWvXnkEC0Btt3MpVhn41q0ddL4qtZ1g9l2ooeVweHH9X5xCAaG3R28lpzivQzj8AZRVSPHs98yseh0L9I10-G0frVqteUkjuQJQNAmFlEPpqegmDhpPo_vaMScx1w1ivQ--mxoZli2uvBI1Pnen-nGIyD3qT_oT-0r3K6bjwk1vI71fgMtGw6ZhWWFuWv48Oq-vg7lN1fKgPEyvcJMcdIISnra0z-1w/file
#'
#' @export
#'
RForm_VAR<-function(TSL,p,W=NULL){
  for(i in 1:p){
    if(i==1){
      aux<-lag.xts(TSL,k=i)
    }else{
      aux<-cbind(aux,lag.xts(TSL,k=i))
    }
  }
  Y<-TSL[(p+1):nrow(TSL),]

  if(!is.null(W)){
    X    = cbind(W[(p+1):nrow(W),],aux[(p+1):nrow(aux),])

    m    = nrow(W)
  }else{
    X    = cbind(matrix(1,nrow(Y),1),aux[(p+1):nrow(aux),])
    m    = 1
  }
  slopeparameters = t(solve(t(X)%*%X)%*%t(X)%*%Y)
  t(solve(crossprod(X), crossprod(X,Y)))

  AL   = slopeparameters[,(m+1):ncol(slopeparameters)]
  mu   = slopeparameters[,1:m]

  eta  = t(Y)-slopeparameters%*%(t(X))

  Sigma= (eta%*%t(eta))/(ncol(eta))

  output<-list(mu=mu,AL=AL,Sigma=Sigma,eta=eta,X=X,Y=Y,p=p)
  return(output)
}
