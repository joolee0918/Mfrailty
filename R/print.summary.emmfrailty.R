# ========================================================= 
# Print class summary.emmfrailty 
# =========================================================

#' @export
#' @method print summary.emmfrailty
#' @keywords internal
print.summary.emmfrailty <- function(x, digits = max(getOption("digits") - 3, 3), 
    signif.stars = getOption("show.signif.stars"), ...) {
    
    if (!is.null(x$fail)) {
        cat("Fitting failed.", x$fail, "\n")
        return()
    }
    
    savedig <- options(digits = digits)
    on.exit(options(savedig))
    J = x$neventtype
    
    K = J * (J - 1)/2
    com = combn(J, 2)
    tmp <- list()
    
    cat("\n")
    if (is.null(x$copula)) 
        cat("#### Semiparametric Independent Estimation ####") else {
        if (x$two_stage) 
            cat("#### Semiparametric Two-stage Pairwise Estimation ####") else cat("#### Semiparametric Pairwise Estimation ####")
    }
    cat("\n\n")
    
    cat("Marginal Parameters")
    for (i in 1:J) {
        cat("\n")
        cat(" Event: ", x$varnames$eventnames[i])
        cat("\n")
        if (!is.null(x$coefficients[[i]])) {
            printCoefmat(x$coefficients[[i]], signif.stars = signif.stars, ...)
        }
        if (!is.null(x$conf.int[[i]])) {
            cat("\n")
            print(x$conf.int[[i]])
        }
        cat("\n")
        cat("Variance of random effect ", x$margins[i], ": ", x$sig2[i], " ", "(", 
            sqrt(x$sig2var[i]), ")", sep = "", "\n")
        cat("n=", x$n)
        if (!is.null(x$nevent[i])) 
            cat(", number of events", x$varnames$eventnames[i], "=", x$nevent[i], 
                "\n") else cat("\n")
        cat("--------------------------------------------------------------")
    }
    if (!is.null(x$copula)) {
        cat("\n\n")
        cat("Depedente Parameter: ", x$copula, sep = "", "\n")
        for (i in 1:K) {
            cat("Event", com[1, i], "&", "Event", com[2, i], "- ")
            if (x$copula == "student") {
                cat("Copula parameter ", ": ", x$rho[i], " ", "(", sqrt(x$rhovar[i]), 
                  ")", sep = "")
                cat(" ", x$nu[i], " ", "(", sqrt(x$nuvar[i]), ")", sep = "")
            } else {
                cat("Copula parameter ", ": ", x$rho[i], " ", "(", sqrt(x$rhovar[i]), 
                  ")", sep = "")
            }
            cat("\n")
        }
    }
    cat("\n")
    cat("Observed log Likelihood:", x$logLik)
    cat("\nIntegration: (Gauss-Legendre Quadrature) ")
    cat("quadrature points:", x$control$nknot)
    if (is.null(x$copula)) {
        cat("\nConvergence: \n")
        for (i in 1:J) cat(" Event", x$varnames$eventnames[i], ": ", x$convergence[i], 
            sep = "")
    } else {
        if (x$two_stage) {
            cat("\nConvergence: \n")
            cat(" First Stage: ")
            for (i in 1:J) cat(" Event", x$varnames$eventnames[i], ": ", x$convergence[i], 
                sep = "")
            cat("\n Second Stage:", x$convergence[J + 1])
        } else {
            cat("\nConvergence:", x$convergence)
        }
    }
    cat("\n")
    invisible(x)
    
}
