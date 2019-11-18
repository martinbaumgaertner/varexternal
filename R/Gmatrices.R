#' Gmatrices
#'
#' Computes the derivatives of vec(C) wrt vec(A) based on Lütkepohl H. New introduction to multiple time series analysis. Springer, 2007.
#'
#' @param AL VAR model coefficients
#' @param C MA representation coefficients
#' @param p lag order
#' @param hori forecast horizon
#' @param n number of variables
#'
#' @return G: derivatives of C wrt A
#'
#' @seealso https://pdfs.semanticscholar.org/3e18/3a5ec97ff636363e4deedff7eaeee9d894c9.pdf
#'
#' @export
Gmatrices<-function(AL,C,p,hori,n){
  J = cbind(eye(n), zeros(n,(p-1)*n))
  Alut = rbind(AL,
               cbind(eye(n*(p-1)),zeros(n*(p-1),n))
               )
  AJ=array(rep(zeros(n*p, n), hori), c(n*p, n, hori))
  for (k in 1:hori){
    AJ[,,k] = ((Alut)%^%(k-1)) %*% t(J)
  }
  JAp = t(array(AJ, c(n*p,n*hori)))
  AJaux = array(rep(zeros(size(JAp,1)*n, size(JAp,2)*n), hori), c(size(JAp,1)*n, size(JAp,2)*n, hori))

  Caux = array(cbind(eye(n), C[,(1:((hori-1)*n))]), c(n,n,hori))

  for (i in 1:(hori)){
    AJaux[(((n^2)*(i-1))+1):dim(AJaux)[1],,i] = kron(JAp[1:(n*(hori+1-i)),], Caux[,,i])
  }

  Gaux = aperm(array(t(apply(AJaux,1:2,sum)), c((n^2)*p, n^2, hori)), c(2,1,3))
  G = array(rep(zeros(dim(Gaux)[1], dim(Gaux)[2])),c(dim(Gaux)[1], dim(Gaux)[2], dim(Gaux)[3]+1))
  G[,,2:dim(G)[3]] = Gaux;

  for (i in 2:dim(G)[3]){
    if(i==2){
      Gcum<-cbind(G[,1:dim(G)[2],2])
    }else{
      Gcum<-cbind(Gcum,G[,1:dim(G)[2],i])
    }
  }
  Gcum <- array(Gcum, c(dim(G)[1], dim(G)[2], dim(G)[3]))

  return(list(G=G,
              Gcum=Gcum))
}




