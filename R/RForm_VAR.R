#' RForm_VAR
#'
#' Provides reduced form estimators of a VAR(p) model
#' @param TSL matrix of time series
#' @param p number of lags in the VAR model
#' @param W Matrix of exogenous regressors
#' @return AL Least-squares estimator of the VAR coefficients
#' @return Sigma Least-squares estimator of the VAR residuals
#' @return eta VAR model residuals
#' @return X Matrix of VAR covariates
#' @return Y VAR matrix of endogenous regressors
#'
#' If no exogenous regressors are specified, our estimation always includes a constant.
#'
#' @seealso https://github.com/jm4474/SVARIV
#' @export
RForm_VAR<-function(TSL,p,W=NULL){
  for(i in 1:p){
    if(i==1){
      aux<-xts::lag.xts(TSL,k=i)
    }else{
      aux<-cbind(aux,xts::lag.xts(TSL,k=i))
    }
  }
  Y<-TSL[(p+1):nrow(TSL),]

  if(!is.null(W)){
    X = cbind(W[(p+1):nrow(W),],aux[(p+1):nrow(aux),])
    m = nrow(W)
  }else{
    X = cbind(matrix(1,nrow(Y),1),aux[(p+1):nrow(aux),])
    m = 1
  }
  slopeparameters = t(solve(t(X)%*%X)%*%t(X)%*%Y)

  AL = slopeparameters[,(m+1):ncol(slopeparameters)]
  mu = slopeparameters[,1:m]

  eta = t(Y)-slopeparameters%*%(t(X))

  Sigma= (eta%*%t(eta))/(ncol(eta))

  output<-list(mu=mu,AL=AL,Sigma=Sigma,eta=eta,X=X,Y=Y,p=p)
  return(output)
}
