#' lee.R
#' Lee. et al's method of approximating the mixture distribution of Q to match the theoretical distribution
#' this is based on Liu. et al's method
#' Helper function for getting the parameters for the test statistics approximation
#' @importFrom bindata bincorr2commonprob commonprob2sigma rmvbin
#' @importFrom stats dbeta rnorm qchisq pchisq rbinom
#' @importFrom SKAT SKAT SKAT_Null_Model
#' @param c A vector of 4 elements c=c(c1,c2,c3,c4).
#'
#' @return A list with the elements:
#' \item{l}{The degree of freedom of chi-square distribution.}
#' \item{delta}{The non-centroality parameter of chi-square distribution.}
#' \item{mux}{The mean of chi-square distribution}
#' \item{sigmax}{The standard deviation of chi-square distribution}
#' \item{muq}{The mean of test statistics.}
#' \item{muq}{The standard deviation of test statistics.}
#' @export
#' @examples
#' c<-c(1,2,3,4)
#' lee(c)
#' # Should return l=2.220446e-15, delta=1, mux=1, sigmax=2, muq=1, sigmaq=2.

lee<-function(c){
  muq=c[1]
  sigmaq=sqrt(2*c[2])
  s1=c[3]/c[2]^(3/2)
  s2=c[4]/c[2]^2
  if (s1^2>s2)
  {a=1/(s1-sqrt(s1^2-s2))
  delta= (s2*a^4-a^2)/2
  }else
  {a=1/sqrt(s2)
  delta=0}
  l=a^2-2*delta
  mux=l+delta
  sigmax=sqrt(2)*sqrt(l+2*delta)
  return(list(l=l,delta=delta,mux=mux,sigmax=sigmax,muq=muq,sigmaq=sigmaq))
}
