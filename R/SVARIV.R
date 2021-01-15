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
#' @param ci_type     confidence intervals to include (choose from "msw" (Montiel, Stock and Watson, 2020), "delta" or "plugin")
#' @param print_wald     Number of horizons for the Impulse Response Functions (does not include the impact or horizon 0)
#' @param instrument_name     Number of horizons for the Impulse Response Functions (does not include the impact or horizon 0)
#'
#'
#' @return irfs:       list containing all irf data
#' @return waldstat: contain msw waldstat
#'
#' @examples
#'
#'
#'
#' p           = 24    #Number of lags in the VAR model
#' NWlags      = 0;  # Newey-West lags(if it is neccessary to account for time series autocorrelation)
#' norm        = 1; # Variable used for normalization
#' scale       = 1; # Scale of the shock
#' horizons    = 20; #Number of horizons for the Impulse Response Functions(IRFs)
#' confidence=c(0.6,0.9,0.95);
#' data(oil)
#' colnames(oil)<-c("a","b","c","d","year","month")
#' ydata<-oil[,1:3]
#' z<-oil[,4]
#' VAR<-SVARIV(ydata,z,p,confidence,NWlags,norm,scale,horizons,instrument_name="test")
#'
#' @export

SVARIV<-function(ydata, z, p, confidence, NWlags, norm, scale, horizons,ci_type=c("msw"),print_wald=T,instrument_name){
  SVARIV_Check(p,confidence, ydata, z, NWlags, norm, scale, horizons)

  SVARinp<-list(ydata=ydata,
                Z=z,
                n=ncol(ydata))

  RForm<-RForm_VAR(SVARinp$ydata, p)
  RForm$Gamma<-RForm$eta%*%SVARinp$Z[(p+1):length(SVARinp$Z)]/ncol(RForm$eta)
  RForm$Y0         = SVARinp$ydata[1:p,]
  RForm$externalIV = SVARinp$Z[(p+1):length(SVARinp$Z)]
  RForm$n          = SVARinp$n
  RForm$names<-colnames(ydata)
  n=RForm$n
  Ti=nrow(RForm$eta)
  d=((n^2)*p)+(n)
  dall= d+ (n*(n+1))/2

  RForm<-c(RForm,CovAhat_Sigmahat_Gamma(p,RForm$X,SVARinp$Z[(p+1):length(SVARinp$Z)],RForm$eta,NWlags))

  InferenceMSW = MSWfunction(confidence,norm,scale,horizons,RForm)

  if(any(ci_type=="msw")){
    waldstat<-melt.list(InferenceMSW) %>%
      filter(L3 %in% c("critval","Waldstat")) %>%
      select(-X1,-X2,-L2) %>%
      dplyr::rename("type"=L3,"level"=L1) %>%
      unique() %>%
      filter(row_number()!=1) %>%
      arrange(value)

    position=max(which(waldstat$type=="Waldstat"))

    if(position!=length(waldstat$type)){
      if (print_wald) {
        message('NOTE: The Wald statistic for the covariance between the instrument and the normalized variable is:')
        message(round(waldstat %>%
                        filter(type=="Waldstat") %>% select(value) %>% pull(),3))
        message('Given the confidence level, if the Wald statistic is larger than: ',round(waldstat %>%
                                                                                             filter(type=="critval") %>% select(value) %>% pull(),3),"(",confidence,")")
        message('the weak-IV robust confidence set will be a bounded interval for every horizon (check "casedummy" if not).')
      }
    }
    }else{
      Stat<-NULL
    }

  irfs=as_tibble(melt.list(InferenceMSW) %>%
    dplyr::rename("variable"=X1,
                  "horizon"=X2,
                  "type"=L3,
                  "confi_type"=L2,
                  "confi_level"=L1) %>%
    filter(complete.cases(.))) %>%
    mutate(instrument=instrument_name)

  return(list(
    irfs=irfs,
    waldstat=waldstat
  ))
}
