#' Point estimates for two-sided M-Wright distribution
#'
#'
#' This provides point estimates for the shape and scale parameters.
#'
#'
#' @param x numeric vector.
#'
#'
#' @return numeric vector
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
#
#' x <- rmwright2(1000, 0.7, 0.4)
#' point_est2(x)
#'
#'
#' @import stats
#'
#' @export
point_est2 <- function(x) {
    ec=0.57721566490153286
    x=x-mean(x)
    x=x[x!=0]
    x <- log(abs(x))
    ah<- ( ( 1- 6*var(x) /pi^2 )^2)^(1/4)
    sh <- exp(ec * (1 - ah) + mean(x))
    return( c(ah, sh) )
}
