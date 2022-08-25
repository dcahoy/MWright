#' Interval estimation  for one-sided M-Wright distribution
#'
#'
#' Confidence intervals for the model parameters.
#'
#'
#' @param x numeric vector
#' @param lev confidence level.
#'
#'
#'
#' @return matrix
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
#' mwright_1sided <- rmwright1(1000, 0.7, 0.4)
#' int_est1(mwright_1sided ,0.95)
#'
#'
#' @import stats
#'
#'
#' @export
int_est1 <- function(x, lev) {
  ec <- 0.57721566490153286
  zeta3 <- 1.202056903159594285399
  dum<-point_est1(x)
  n<-length(x)
  ah<-dum[1]
  sh<-dum[2]
  zcv <- qnorm(1 - (1 - lev)/2, 0, 1)
  sah <- sqrt((-1 - (ah^4 - 11)/(10 * ah^2))/n)
  ssh <- sqrt((((sh^2) * (-(-1 + ah^2) * (pi^2) * (3 * (11 + ah^2) *ec^2 + 5 * (ah^2) * (pi^2))
                         + 360 * ah * (-1 + ah^3) * ec *(zeta3)))/(30 * (ah^2) * (pi^2)))/n)
  return( rbind( c('alpha:', ah - zcv * sah, ah + zcv * sah), c('rho:', sh - zcv * ssh, sh + zcv * ssh) ) )
}



