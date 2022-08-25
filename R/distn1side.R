#' Distribution function for one-sided M-Wright distribution
#'
#'
#' Calculates a left-tail probability.
#'
#'
#' @param alp point estimate for shape parameter alpha.
#' @param sc point estimate for scale parameter s.
#' @param upper non-negative upper quantile
#'
#'
#' @return numeric
#'
#'
#' @references
#'
#' Cahoy and Minkabo (2017). \emph{Inference for three-parameter M-Wright distributions with applications.} Model Assisted Statistics and Applications, 12(2), 115-125.
#' \url{https://doi.org/10.3233/MAS-170388}
#'
#' Cahoy (2012). \emph{Moment estimators for the two-parameter M-Wright distribution.} Computational Statistics, 27(3), 487-497.
#' \url{https://doi.org/10.1007/s00180-011-0269-x}
#'
#' Cahoy (2012). \emph{Estimation and simulation for the M-Wright function.}  Communications in Statistics-Theory and Methods, 41(8), 1466-1477.
#' \url{https://doi.org/10.1080/03610926.2010.543299}
#'
#' Cahoy (2011). \emph{On the parameterization of the M-Wright function.} Far East Journal of Theoretical Statistics, 34(2), 155-164.
#' \url{http://www.pphmj.com/abstract/5767.htm}
#'
#' Mainardi, Mura, and Pagnini (2010). \emph{The M-Wright Function in Time-Fractional Diffusion Processes: A Tutorial Survey}. Int. J. Differ. Equ., Volume 2010.
#' \url{https://doi.org/10.1155/2010/104505}
#'
#'
#' @examples
#'
#' pmwright1(runif(1), runif(1,0,10),Inf  )
#'
#' pmwright1(runif(1), runif(1,0,10), 0.5 )
#'
#'
#'
#' @import stats cubature
#'
#'
#' @export
pmwright1 <-function(alp,sc, upper){
  b<-alp
  s<-sc
  f<- function(x){
    (1/s)*exp(-(((x[2]/s)^(-1/b))^(b/(b-1)))*((sin(b*(x[1]+pi/2))/cos(x[1]))^(b/(1-b)))*cos(x[1]*(b-1)+b*pi/2)/cos(x[1]))*((sin(b*(x[1]+pi/2))/cos(x[1]))^(b/(1-b)))*(cos(x[1]*(b-1)+b*pi/2)/cos(x[1]))*(b*((x[2]/s)^(-1/b))^(1/(b-1)))/(pi*abs(1-b))*( (1/b)*( (x[2]/s)^(-1-1/b)) )
  }
  return( adaptIntegrate(f, lowerLimit = c(-pi/2, 0), upperLimit = c(pi/2, upper))$integral )
}



