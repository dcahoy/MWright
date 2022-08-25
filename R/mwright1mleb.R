#' One-sided M-Wright distribution
#'
#'
#' Plots the density function.
#'
#'
#' @param ah point estimate for shape parameter alpha.
#' @param sh point estimate for scale parameter s.
#' @param  m number of data points (pairs) to use for plotting.
#' @param  max maximum  x-axis value to use for plotting.
#'
#'
#' @return numeric matrix
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
#' xy=dmwright1(0.45, 2.5, 1000, 10)
#' plot(xy[,1], xy[,2],  lwd = 2, type="l",ylab="", xlab="x")
#'
#' mwright1_sided <- rmwright1(1000, 0.45, 2.5)
#' hist(mwright1_sided, br=30, prob=TRUE)
#' lines(xy[,1], xy[,2],  lwd=2 )
#'
#'
#'
#' @import stats
#'
#'
#' @export
dmwright1b<- function(x, para[1] ,para[2]) {
      s<-para[2]
      b <- para[1]
      y<-(x/s)^(-1/b)
      int <- function(u){
        exp(-(y^(b/(b-1)))*((sin(b*(u+pi/2))/cos(u))^(b/(1-b)))*cos(u*(b-1)+b*pi/2)/cos(u))*((sin(b*(u+pi/2))/cos(u))^(b/(1-b)))*cos(u*(b-1)+b*pi/2)/cos(u)
      }
      smot<-(1/s)*(integrate(int, lower = -pi/2, upper = pi/2)$value)*(b*y^(1/(b-1)))/(pi*abs(1-b))*( (1/b)*( (x/s)^(-1-1/b)) )

    return( smot )
}



