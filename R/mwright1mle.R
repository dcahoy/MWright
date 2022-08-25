#' Maximum likelihood estimates for one-sided M-Wright distribution
#'
#'
#' Maximum likelihood estimates
#'
#'
#' @param n  sample size.
#' @param nu a number between 0 and 1.
#' @param sc a non-negative scale value.
#'
#'
#'
#' @return MLEs for shape and scale parameters
#'
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
#'
#' @examples
#'
#'
#'  ????
#'
#'
#' @import stats
#'
#'
#' @export
mw1mle <-function

dmwright1b<- function(x, par) {
  s<- par[2]
  b <- par[1]
  y<-(x/s)^(-1/b)
  int <- function(u){
    exp(-(y^(b/(b-1)))*((sin(b*(u+pi/2))/cos(u))^(b/(1-b)))*cos(u*(b-1)+b*pi/2)/cos(u))*((sin(b*(u+pi/2))/cos(u))^(b/(1-b)))*cos(u*(b-1)+b*pi/2)/cos(u)
  }
   (1/s)*(integrate(int, lower = -pi/2, upper = pi/2)$value)*(b*y^(1/(b-1)))/(pi*abs(1-b))*( (1/b)*( (x/s)^(-1-1/b)) )

}

dmwright1b<- function(x, par) {
  s<- par[2]
  b <- par[1]
  y<-(x/s)^(-1/b)
  int <- function(u){
    exp(-(y^(b/(b-1)))*((sin(b*(u+pi/2))/cos(u))^(b/(1-b)))*cos(u*(b-1)+b*pi/2)/cos(u))*((sin(b*(u+pi/2))/cos(u))^(b/(1-b)))*cos(u*(b-1)+b*pi/2)/cos(u)
  }
  (1/s)*(integrate(int, lower = -pi/2, upper = pi/2)$value)*(b*y^(1/(b-1)))/(pi*abs(1-b))*( (1/b)*( (x/s)^(-1-1/b)) )

}
llike <-function(x,para){
return(-sum(log(dmwright1b(x,para))))
}
x <- rmwright1(1000, 0.7, 0.4)

optim(par=point_est1(x), llike, x=x)


