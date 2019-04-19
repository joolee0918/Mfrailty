# ========================================================= 
# Control List
# =========================================================

#' Control parameters for emmfrailty and parmfrailty

#' @param tol.P Tolerance for convergence in the parameters. Default is 1e-04.
#' @param tol.L Tolerance for convergence in the log-likelihood. Default is 1e-10.
#' @param maxiter The maxiumum number of interations in Expectation-Maximization Algorithm for a semiparametric model and in \code{optim} for a parametric model.
#' @param delta maxiter An increament used for numerical differentiation in profile likelihood. Default is 1e-03.
#' @param nknot The number of nodes used in Gaussian-quadrature. Default is 20.

controlList <- function(control) {
    
    control.default <- list(tol.P = 1e-04, tol.L = 1e-10, maxiter = 2000, delta = 1e-03, 
        nknot = 20)
    
    control <- c(control)
    namec <- names(control)
    
    if (length(uname <- namec[!namec %in% names(control.default)]) > 0) {
        warning("\n unknown names in 'control': ", paste(uname, collapse = ", "))
    }
    
    control.default[namec] <- control
    
    return(control.default)
    
}




