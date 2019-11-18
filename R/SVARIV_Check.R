#' SVARIV_Check
#'
#' Checks whether the inputs from SVARIV are valid.
#'
#' @param p Number of lags in the VAR model
#' @param confidence   Value for the standard and weak-IV robust confidence set
#' @param ydata        Endogenous variables from the VAR model
#' @param z            External instrumental variable
#' @param NWlags       Newey-West lags
#' @param norm         Variable used for normalization
#' @param scale        Scale of the shock
#' @param horizons     Number of horizons for the Impulse Response Functions
#'
#' @return None
#'
#'
#' @export
SVARIV_Check<-function(p,confidence, ydata, z, NWlags, norm, scale, horizons){
  Ti=nrow(ydata)
  if(is.null(p)){
    stop('p must be assigned a value.')
  }
  if(!is.numeric(p)){
    stop('p must be numeric.')
  }
  if(length(p)!=1){
    stop('p must have only one element')
  }
  if(is.null(confidence)){
    stop('confidence must be assigned a value.')
  }
  if(!is.numeric(confidence)){
    stop('confidence must be numeric.')
  }
  if(length(confidence)!=1){
    stop('confidence must have only one element')
  }
  if(confidence>=1|confidence<=0){
    stop('confidence must be > 0 and < 1')
  }
  if(is.null(ydata)){
    stop('ydata must be a (T times n) matrix. It is now empty')
  }
  if(!is.matrix(ydata)){
    stop('ydata must be a matrix')
  }
  if(nrow(ydata)<ncol(ydata)){
    stop('ydata: number of rows < number of columns. ydata must be T times n.')
  }
  if(is.null(z)){
    stop('z must be (T times 1) vector. It is now empty')
  }
  if(ncol(z)!=1){
    stop('There must be only one instrument (z must have only one column)')
  }
  if(nrow(z)!=T){
    stop('z must have T rows')
  }
}
