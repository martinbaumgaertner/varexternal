#' SvAR-IV.
#'
#' Implements standard and weak-IV robust SVAR-IV inference.
#'
#' @param p Number of lags in the VAR model
#' @param confidence   Value for the standard and weak-IV robust confidence set
#' @param ydata        Endogenous variables from the VAR model
#' @param z            External instrumental variable
#' @param NWlags       Newey-West lags
#' @param norm         Variable used for normalization
#' @param scale        Scale of the shock
#' @param horizons     Number of horizons for the Impulse Response Functions
#' @param savdir       Directory where the figures generated will be saved
#' @param columnnames  Vector with the names for the endogenous variables, in the same order as ydata
#' @param IRFselect    Indices for the variables that the user wants separate IRF plots for
#' @param cumselect    Indices for the variables that the user wants cumulative IRF plots for
#' @param time         Time unit for the dataset e.g. year, month, etc.
#' @param dataset_name The name of the dataset used for generating the figures
#'
#' @return      PLugin:       Structure containing standard plug-in inference
#' @return       InferenceMSW: Structure containing the MSW weak-iv robust confidence interval
#' @return       Chol:         Cholesky IRFs
#' @return       RForm:        Structure containing the reduced form parameters
#'
#' @export

SVARIV<-function(p,confidence, ydata, z, NWlags, norm, scale, horizons, savdir, columnnames, IRFselect,
                 cumselect, time, dataset_name){
  SVARIV_Check(p,confidence, ydata, z, NWlags, norm, scale, horizons, savdir, columnnames, IRFselect, cumselect, time, dataset_name)

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

