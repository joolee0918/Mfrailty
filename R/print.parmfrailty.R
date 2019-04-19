# ========================================================= 
# Print class parmfrailty 
# =========================================================

#' @export
#' @method print parmfrailty
#' @keywords internal
print.parmfrailty <- function(object, ...) {
    
    if (!is.null(object$fail)) {
        cat("Fitting failed.", object$fail, "\n")
        return()
    }
    
    digits = max(1L, getOption("digits") - 3L)
    signif.stars = FALSE
    savedig <- options(digits = digits)
    on.exit(options(savedig))
    J <- object$neventtype
    
    beta.coef <- object$theta.coefficients
    se <- sqrt(diag(object$var))
    beta.se <- lapply(1:J, function(i) sqrt(object$thetavar[[i]]))
    if (is.null(beta.coef) | is.null(beta.se)) 
        stop("Input is not valid")
    
    K = J * (J - 1)/2
    com = combn(J, 2)
    tmp <- list()
    for (i in 1:J) {
        if (!is.null(object$betavar[[i]])) {
            tmp[[i]] <- cbind(beta.coef[[i]], exp(beta.coef[[i]]), beta.se[[i]], 
                beta.coef[[i]]/beta.se[[i]], pchisq((beta.coef[[i]]/beta.se[[i]])^2, 
                  1, lower.tail = FALSE))
            dimnames(tmp[[i]]) <- list(names(beta.coef[[i]]), c("coef", "exp(coef)", 
                "se(coef)", "z", "p"))
        }
    }
    
    cat("\n")
    cat("#### Parametric Pairwise Estimation ####")
    
    cat("\n\n")
    cat("Marginal Parameters")
    for (i in 1:J) {
        cat("\n")
        cat("Event: ", object$varnames$eventnames[i], object$dist[i])
        cat("\n")
        printCoefmat(tmp[[i]], digits = digits, P.values = TRUE, has.Pvalue = TRUE, 
            signif.stars = signif.stars, ...)
        cat("\n")
        cat("Variance of random effect ", object$margins[i], ": ", object$sig2[i], 
            " ", "(", sqrt(object$sig2var[i]), ")", sep = "")
        cat("  n=", object$n)
        if (!is.null(object$nevent[i])) 
            cat(", number of events", object$varnames$eventnames[i], "=", object$nevent[i], 
                "\n") else cat("\n")
    }
    cat("\n")
    cat("\n")
    cat("Depedente Parameter: ", object$copula, sep = "", "\n")
    for (i in 1:K) {
        cat("Event", com[1, i], "&", "Event", com[2, i], "- ")
        cat("Copula parameter ", ": ", object$rho[i], " ", "(", sqrt(object$rhovar[i]), 
            ")", sep = "")
        cat("\n")
    }
    
    invisible(object)
}
