#' Random number generation for two-sided M-Wright distribution
#'
#'
#' Generates random numbers.
#'
#'
#' @param n sample size.
#' @param nu a number between 0 and 1.
#' @param sc a non-negative scale value.
#'
#'
#'
#' @return a vector of two-sided M-Wright distributed random numbers.
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
#' Cahoyd (2011). \emph{On the parameterization of the M-Wright function.} Far East Journal of Theoretical Statistics, 34(2), 155-164.
#' \url{http://www.pphmj.com/abstract/5767.htm}
#'
#' Mainardi, Mura, and Pagnini (2010). \emph{The M-Wright Function in Time-Fractional Diffusion Processes: A Tutorial Survey}. Int. J. Differ. Equ., Volume 2010.
#' \url{https://doi.org/10.1155/2010/104505}
#'
#'
#'
#' @examples
#'
#' mwright_2sided <- rmwright2(1000, 0.7, 0.4)
#' hist(mwright_2sided, br=30)
#'
#' @import stats
#'
#' @export
rmwright2 <- function(n, nu, sc){
    U.1 <- runif(n);
    U.2 <- runif(n);
    U.3 <- runif(n);
    rstab <- ( ( sin(nu*pi*U.2)*sin((1-nu)*pi*U.2)^(1/nu-1) )/ (  ( sin(pi*U.2)^(1/nu) )*abs(log(U.3))^(1/nu-1))  );
    return(sc*( rstab^(-nu) )*ifelse(runif(n)<=0.5, -1, 1));
}


