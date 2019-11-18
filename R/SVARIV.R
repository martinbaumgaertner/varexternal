#' SvAR-IV.
#'
#' Implements standard and weak-IV robust SVAR-IV inference.
#'
#' @param ydata        Endogenous variables from the VAR model
#' @param z            External instrumental variable
#' @param p Number of lags in the VAR model
#' @param confidence   Value for the standard and weak-IV robust confidence set
#' @param NWlags       Newey-West lags (set it to 0 to compute heteroskedasticity robust std errors)
#' @param norm         Variable used for normalization
#' @param scale        Scale of the shock
#' @param horizons     Number of horizons for the Impulse Response Functions (does not include the impact or horizon 0)
#'
#' @return      PLugin:       Structure containing standard plug-in inference
#' @return       InferenceMSW: Structure containing the MSW weak-iv robust confidence interval
#' @return       Chol:         Cholesky IRFs
#' @return       RForm:        Structure containing the reduced form parameters
#'
#' @examples
#' p           = 24    #Number of lags in the VAR model
#' confidence  = .95    #Confidence Level for the standard and weak-IV robust confidence set
#' NWlags      = 0;  # Newey-West lags(if it is neccessary to account for time series autocorrelation)
#' norm        = 1; # Variable used for normalization
#' scale       = 1; # Scale of the shock
#' horizons    = 20; #Number of horizons for the Impulse Response Functions(IRFs)
#' data(oil)
#' ydata<-data[,3:5]
#' z<-data[,6]
#' d<-SVARIV(ydata, z, p, confidence,  NWlags, norm, scale, horizons)
#'
#' @export

SVARIV<-function(ydata, z, p, confidence, NWlags, norm, scale, horizons){
  SVARIV_Check(p,confidence, ydata, z, NWlags, norm, scale, horizons)

  SVARinp<-list(ydata=ydata,
                Z=z,
                n=ncol(ydata))

  RForm<-RForm_VAR(SVARinp$ydata, p)
  RForm$Gamma<-RForm$eta%*%SVARinp$Z[(p+1):nrow(SVARinp$Z),1]/ncol(RForm$eta)
  RForm$Y0         = SVARinp$ydata[1:p,]
  RForm$externalIV = SVARinp$Z[(p+1):nrow(SVARinp$Z),1]
  RForm$n          = SVARinp$n

  n=RForm$n
  Ti=nrow(RForm$eta)
  d=((n^2)*p)+(n)
  dall= d+ (n*(n+1))/2

  RForm<-c(RForm,CovAhat_Sigmahat_Gamma(p,RForm$X,SVARinp$Z[(p+1):nrow(SVARinp$Z),1],RForm$eta,NWlags))

  InferenceMSW = MSWfunction(confidence,norm,scale,horizons,RForm,1)$InferenceMSW
  Plugin = MSWfunction(confidence,norm,scale,horizons,RForm,0)$Plugin
  Chol = MSWfunction(confidence,norm,scale,horizons,RForm,0)$Chol

  return(list(
    InferenceMSW=InferenceMSW,
    Plugin=Plugin,
    Chol=Chol
  ))
}

