#' MARep
#'
#'Reports the confidence interval for IRF coefficients described in Montiel-Olea, Stock, and Watson (2017).
#'
#' @param confidence confidence level
#' @param nvar variable defining the normalization
#' @param scale scale of the shock
#' @param horizons Number of horizons for the IRF
#' @param RForm reduced-form structure
#' @param display_on dummy variable
#'
#' @return InferenceMSW: Structure containing the MSW weak-iv robust confidence interval
#' @return Plugin: Structure containing standard plug-in inference
#' @return Chol: Cholesky IRFs
#'
#' @seealso https://github.com/jm4474/SVARIV
#' @seealso https://uc4384f22718973d3d912d614f2e.dl.dropboxusercontent.com/cd/0/inline2/Aslwi_FA9z6v37dfxMeUP3mEaP-a4ZA8RFrXd5h_EymM1lqGy-UUwGJB3TxqRxluMhE9ebnvHEMZktfcUmJ1v4mGSxeS_ZEDJBCrjGvyopAnP44ClPLJ6kmM36F557kHnqKQjn___0XTNvPaWvXnkEC0Btt3MpVhn41q0ddL4qtZ1g9l2ooeVweHH9X5xCAaG3R28lpzivQzj8AZRVSPHs98yseh0L9I10-G0frVqteUkjuQJQNAmFlEPpqegmDhpPo_vaMScx1w1ivQ--mxoZli2uvBI1Pnen-nGIyD3qT_oT-0r3K6bjwk1vI71fgMtGw6ZhWWFuWv48Oq-vg7lN1fKgPEyvcJMcdIISnra0z-1w/file
#'
#' @export
#'
MSWfunction <-function(confidence,nvar,scale,horizons,RForm,display_on) {
    critval = qnorm(1 - ((1 - confidence) / 2), 0, 1) ^ 2
    Caux        = cbind(eye(RForm$n), MARep(RForm$AL, RForm$p, horizons))
    C = array(Caux, c(RForm$n, RForm$n, horizons + 1))
    for (i in 2:dim(C)[3]) {
      if (i == 2) {
        Ccum <- cbind(C[, 1:dim(C)[2], 2])
      } else{
        Ccum <- cbind(Ccum, C[, 1:dim(C)[2], i])
      }
    }
    Ccum <- array(Ccum, c(dim(C)[1], dim(C)[2], dim(C)[3]))

    G <-Gmatrices(RForm$AL,
                MARep(RForm$AL, RForm$p, horizons),
                RForm$p,
                horizons,
                RForm$n)$G
    Gcum <-Gmatrices(RForm$AL,
                MARep(RForm$AL, RForm$p, horizons),
                RForm$p,
                horizons,
                RForm$n)$Gcum

    B1chol      = t(chol(RForm$Sigma))
    B1chol      = scale * (B1chol[, 1] / B1chol[nvar, 1])

    Chol <- array(numeric(), c(0, 0, 2))
    Chol[, , 1] = colSums(aperm(sweep(C, 2, t(B1chol), FUN = "*"), c(2, 1, 3)))

    Chol[, , 2] = colSums(aperm(sweep(Ccum, 2, t(B1chol), FUN = "*"), c(2, 1, 3)))
    W1          = RForm$WHat[1:((RForm$n ^ 2) * RForm$p), 1:((RForm$n ^ 2) *RForm$p)]
    W12         = RForm$WHat[1:((RForm$n ^ 2) * RForm$p), (1 + (RForm$n ^2) * RForm$p):(ncol(RForm$WHat))]
    W2          = RForm$WHat[(1 + (RForm$n ^ 2) * RForm$p):nrow(RForm$WHat), (1 +(RForm$n ^ 2) * RForm$p):ncol(RForm$WHat)]

    n = RForm$n

    Ti = dim(RForm$eta)[2]

    e = eye(n)
    ahat      = zeros(n, horizons + 1)
    bhat      = zeros(n, horizons + 1)
    chat      = zeros(n, horizons + 1)
    Deltahat  = zeros(n, horizons + 1)
    MSWlbound = zeros(n, horizons + 1)
    MSWubound = zeros(n, horizons + 1)
    casedummy = zeros(n, horizons + 1)
    j=1;ih=2
    for (j in 1:n) {
      for (ih in 1:(horizons + 1)) {
        ahat[j, ih] = (Ti %*% (RForm$Gamma[nvar, 1] ^ 2)) - (critval %*% W2[nvar, nvar])
        bhat[j, ih] = -2 * Ti * scale * (t(e[, j]) %*% C[, , ih] %*% RForm$Gamma) %*%RForm$Gamma[nvar, 1] +
          2 * critval * scale * kronecker(t(RForm$Gamma), t(e[, j])) %*% G[, , ih] %*%W12[, nvar] +
          2 * critval * scale * t(e[, j]) %*% C[, , ih] %*% W2[, nvar]
        chat[j, ih] = ((Ti ^ .5) * scale * t(e[, j]) %*% C[, , ih] %*% RForm$Gamma) ^2 -
          critval * (scale ^ 2) * (kronecker(t(RForm$Gamma), t(e[, j]))) %*%G[, , ih] %*% W1 %*% t((kronecker(t(RForm$Gamma), t(e[, j]))) %*% G[, , ih]) -
          2 * critval * (scale ^ 2) * (kronecker(t(RForm$Gamma), t(e[, j]))) %*%G[, , ih] %*% W12 %*% t(C[, , ih]) %*% e[, j] -
          critval * (scale ^ 2) * t(e[, j]) %*% C[, , ih] %*% W2 %*% t(C[, , ih]) %*%e[, j]

        Deltahat[j, ih] = bhat[j, ih] ^ 2 - (4 * ahat[j, ih] * chat[j, ih])

        if (ahat[j, ih] > 0 && Deltahat[j, ih] > 0) {
          casedummy[j, ih] = 1
          MSWlbound[j, ih] = (-bhat[j, ih] - (Deltahat[j, ih] ^ 0.5)) / (2 *ahat[j, ih])
          MSWubound[j, ih] = (-bhat[j, ih] + (Deltahat[j, ih] ^ 0.5)) / (2 *ahat[j, ih])
        } else if (ahat[j, ih] < 0 && Deltahat[j, ih] > 0) {
          casedummy[j, ih] = 2
          MSWlbound[j, ih] = (-bhat[j, ih] + (Deltahat[j, ih] ^ 0.5)) / (2 *ahat[j, ih])
          MSWubound[j, ih] = (-bhat[j, ih] - (Deltahat[j, ih] ^ 0.5)) / (2 *ahat[j, ih])
        } else if (ahat[j, ih] > 0 && Deltahat[j, ih] < 0) {
          casedummy[j, ih] = 3
          MSWlbound[j, ih] = NA
          MSWubound[j, ih] = NA
        } else{
          casedummy[j, ih] = 4
          MSWlbound[j, ih] = -Inf
          MSWubound[j, ih] = Inf
        }

      }
    }

    MSWlbound[nvar, 1] = scale
    MSWubound[nvar, 1] = scale

    InferenceMSW <- list(
      ahat = ahat,
      bhat = bhat,
      chat = chat,
      Deltahat = Deltahat,
      casedummy = casedummy,
      MSWlbound = MSWlbound,
      MSWubound = MSWubound,
      Ti = Ti
    )

    lambdahat                  = zeros(n, horizons + 1)
    DmethodVar                 = zeros(n, horizons + 1)
    Dmethodlbound              = zeros(n, horizons + 1)
    Dmethodubound              = zeros(n, horizons + 1)

    for (ih in 1:(horizons + 1)) {
      for (ivar in 1:n) {
        lambdahat[ivar, ih]     = scale * t(e[, ivar]) %*% C[, , ih] %*% RForm$Gamma /RForm$Gamma[nvar, 1]
        d1                     = ((kronecker(t(RForm$Gamma), t(e[, ivar])) *scale) %*% G[, , ih])
        d2                     = (scale * t(e[, ivar]) %*% C[, , ih]) -(lambdahat[ivar, ih] %*% t(e[, nvar]))
        d                      = t(cbind(d1, d2))

        DmethodVar[ivar, ih]    = t(d) %*% RForm$WHat %*% d
        Dmethodlbound[ivar, ih] = lambdahat[ivar, ih] - ((critval / Ti) ^0.5) * (DmethodVar[ivar, ih] ^ 0.5) / abs(RForm$Gamma[nvar, 1])
        Dmethodubound[ivar, ih] = lambdahat[ivar, ih] + ((critval / Ti) ^0.5) * (DmethodVar[ivar, ih] ^ 0.5) / abs(RForm$Gamma[nvar, 1])

        rm(d1)
        rm(d2)
        rm(d)

      }
    }

    InferenceMSW <- c(InferenceMSW,
                      list(Dmethodlbound = Dmethodlbound,
                           Dmethodubound = Dmethodubound))

    Plugin <- list(IRF = lambdahat,
                   IRFstderror = (DmethodVar ^ 0.5) / ((Ti ^ 0.5) * abs(RForm$Gamma[nvar, 1])))

    ahatcum      = zeros(n, horizons + 1)
    bhatcum      = zeros(n, horizons + 1)
    chatcum      = zeros(n, horizons + 1)
    Deltahatcum  = zeros(n, horizons + 1)
    MSWlboundcum = zeros(n, horizons + 1)
    MSWuboundcum = zeros(n, horizons + 1)
    casedummycum = zeros(n, horizons + 1)

    for (j in 1:n) {
      for (ih in 1:(horizons + 1)) {
        ahatcum[j, ih] = (Ti %*% (RForm$Gamma[nvar, 1] ^ 2)) - (critval %*% W2[nvar, nvar])
        bhatcum[j, ih] = -2 * Ti * scale * (t(e[, j]) %*% Ccum[, , ih] %*% RForm$Gamma) %*%RForm$Gamma[nvar, 1] +
          2 * critval * scale * kronecker(t(RForm$Gamma), t(e[, j])) %*% Gcum[, , ih] %*%W12[, nvar] +
          2 * critval * scale * t(e[, j]) %*% Ccum[, , ih] %*% W2[, nvar]
        chatcum[j, ih] = ((Ti ^ 0.5) * scale * t(e[, j]) %*% Ccum[, , ih] %*%RForm$Gamma) ^ 2 -
          critval * (scale ^ 2) * (kronecker(t(RForm$Gamma), t(e[, j]))) %*%Gcum[, , ih] %*% W1 %*% t((kronecker(t(RForm$Gamma), t(e[, j]))) %*% Gcum[, , ih]) -
          2 * critval * (scale ^ 2) * (kronecker(t(RForm$Gamma), t(e[, j]))) %*%Gcum[, , ih] %*% W12 %*% t(Ccum[, , ih]) %*% e[, j] -
          critval * (scale ^ 2) * t(e[, j]) %*% Ccum[, , ih] %*% W2 %*% t(Ccum[, , ih]) %*%e[, j]

        Deltahatcum[j, ih] = bhatcum[j, ih] ^ 2 - (4 * ahatcum[j, ih] * chatcum[j, ih])

        if (ahatcum[j, ih] > 0 & Deltahatcum[j, ih] > 0) {
          casedummycum[j, ih] = 1
          MSWlboundcum[j, ih] = (-bhatcum[j, ih] - (Deltahatcum[j, ih] ^ 0.5)) /(2 * ahatcum[j, ih])
          MSWuboundcum[j, ih] = (-bhatcum[j, ih] + (Deltahatcum[j, ih] ^ 0.5)) /(2 * ahatcum[j, ih])
        } else if (ahatcum[j, ih] < 0 & Deltahatcum[j, ih] > 0) {
          casedummycum[j, ih] = 2
          MSWlboundcum[j, ih] = (-bhatcum[j, ih] + (Deltahatcum[j, ih] ^ 0.5)) /(2 * ahatcum[j, ih])
          MSWuboundcum[j, ih] = (-bhatcum[j, ih] - (Deltahatcum[j, ih] ^ 0.5)) /(2 * ahatcum[j, ih])
        } else if (ahatcum[j, ih] > 0 & Deltahatcum[j, ih] < 0) {
          casedummycum[j, ih] = 3
          MSWlboundcum[j, ih] = NA
          MSWuboundcum[j, ih] = NA
        } else{
          casedummycum[j, ih] = 4
          MSWlboundcum[j, ih] = -Inf
          MSWuboundcum[j, ih] = Inf
        }

      }
    }

    MSWlboundcum[nvar, 1] = scale
    MSWuboundcum[nvar, 1] = scale

    InferenceMSW <- c(InferenceMSW,list(ahatcum = ahatcum,
                                        bhatcum = bhatcum,
                                        chatcum = chatcum,
                                        Deltahatcum = Deltahatcum,
                                        casedummycum = casedummycum,
                                        MSWlboundcum = MSWlboundcum,
                                        MSWuboundcum = MSWuboundcum
                                      )
    )

    lambdahatcum                  = zeros(n, horizons + 1)
    DmethodVarcum                 = zeros(n, horizons + 1)
    Dmethodlboundcum              = zeros(n, horizons + 1)
    Dmethoduboundcum              = zeros(n, horizons + 1)

    for (ih in 1:(horizons + 1)) {
      for (ivar in 1:n) {
        lambdahatcum[ivar, ih]     = scale * t(e[, ivar]) %*% Ccum[, , ih] %*% RForm$Gamma /RForm$Gamma[nvar, 1]
        d1                     = ((kronecker(t(RForm$Gamma), t(e[, ivar])) *scale) %*% Gcum[, , ih])
        d2                     = (scale * t(e[, ivar]) %*% Ccum[, , ih]) -(lambdahatcum[ivar, ih] %*% t(e[, nvar]))
        d                      = t(cbind(d1, d2))

        DmethodVarcum[ivar, ih]    = t(d) %*% RForm$WHat %*% d
        Dmethodlboundcum[ivar, ih] = lambdahatcum[ivar, ih] - ((critval / Ti) ^0.5) * (DmethodVarcum[ivar, ih] ^ 0.5) / abs(RForm$Gamma[nvar, 1])
        Dmethoduboundcum[ivar, ih] = lambdahatcum[ivar, ih] + ((critval / Ti) ^0.5) * (DmethodVarcum[ivar, ih] ^ 0.5) / abs(RForm$Gamma[nvar, 1])

        rm(d1)
        rm(d2)
        rm(d)
      }
    }

    InferenceMSW <-c(InferenceMSW,list(Dmethodlboundcum = Dmethodlboundcum,
                                       Dmethoduboundcum = Dmethoduboundcum))

    Plugin <- c(Plugin,list(IRFcum = lambdahatcum,
                   IRFstderrorcum = (DmethodVarcum ^ 0.5) / ((Ti ^ 0.5) * abs(RForm$Gamma[nvar, 1]))))

    Waldstat              = (((Ti ^ 0.5) * RForm$Gamma[nvar, 1]) ^ 2) / RForm$WHat[((n ^
                                                                                       2) * RForm$p) + nvar, ((n ^ 2) * RForm$p) + nvar]

    InferenceMSW <- c(InferenceMSW,
                       list(Waldstat = Waldstat))
    if (display_on == 1) {
      print(
        'NOTE: The Wald statistic for the covariance between the instrument and the normalized variable is:'
      )
      print(Waldstat)
      print('Given the confidence level, if the Wald statistic is larger than: \n\n')
      print(critval)
      print(
        'The weak-IV robust confidence set will be a bounded interval for every horizon (check "casedummy" if not).'
      )
      disp('--')

      #display('Also, the Wald statistic for the covariance between the instrument and the full vector of reduced-form residuals is')
      #WaldstatFull= (Ti)*(RForm.Gamma'*(RForm.WHat(((n^2)*RForm.p)+nvar:end,((n^2)*RForm.p)+nvar:end))^(-1)*RForm.Gamma);
      #display(WaldstatFull);
      #display('The 1-alpha quantile of this statistic is:')
      #display(chi2inv(confidence,n))
    }

    #Plugin<-c(Plugin,
    #          epsilonhat=scale*t(RForm$Gamma)%*%(solve(RForm$Sigma))%*%RForm$eta/RForm$Gamma[nvar,1],
    #          epsilonhatstd = (Plugin$epsilonhat-mean(Plugin$epsilonhat))/std(Plugin$epsilonhat))

    return(list(
      InferenceMSW = InferenceMSW,
      Plugin = Plugin,
      Chol = Chol
    ))
  }
