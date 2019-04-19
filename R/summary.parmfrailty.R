# ========================================================= 
# class summary.parmfrailty 
# =========================================================

#' @export
#' @method summary parmfrailty

summary.parmfrailty <- function(object, conf.int = 0.95, scale = 1, ...) {
    
    if (!is.null(object$fail)) {
        class(object) <- "summary.parmfrailty"
        return(object)
    }
    
    J <- object$neventtype
    beta <- lapply(1:J, function(i) object$theta.coefficients[[i]] * scale)
    if (is.null(object$beta.coefficients)) {
        # Null model
        return(object)  #The summary method is the same as print in this case
    }
    nabeta <- !(is.na(beta))  #non-missing coefs
    beta2 <- beta[nabeta]
    
    
    if (is.null(object$theta) | is.null(object$var)) 
        stop("Input is not valid")
    
    beta.se <- lapply(1:J, function(i) sqrt(object$thetavar[[i]]) * scale)
    
    rval <- list(neventtype = object$neventtype, n = object$n, logLik = object$logLik, 
        margins = object$margins, dist = object$dist, copula = object$copula, iter = object$iter, 
        convergence = object$convergence, control = object$control)
    if (!is.null(object$nevent)) 
        rval$nevent <- object$nevent
    
    K = J * (J - 1)/2
    com = combn(J, 2)
    tmp <- list()
    for (i in 1:J) {
        tmp[[i]] <- cbind(beta[[i]], exp(beta[[i]]), beta.se[[i]], beta[[i]]/beta.se[[i]], 
            pchisq((beta[[i]]/beta.se[[i]])^2, 1, lower.tail = FALSE))
        dimnames(tmp[[i]]) <- list(names(beta[[i]]), c("coef", "exp(coef)", "se(coef)", 
            "z", "Pr(>|z|)"))
    }
    
    rval$coefficients <- tmp
    
    if (conf.int) {
        z <- qnorm((1 + conf.int)/2, 0, 1)
        for (i in 1:J) {
            tmp[[i]] <- cbind(exp(beta[[i]]), exp(-beta[[i]]), exp(beta[[i]] - z * 
                beta.se[[i]]), exp(beta[[i]] + z * beta.se[[i]]))
            dimnames(tmp[[i]]) <- list(names(beta[[i]]), c("exp(coef)", "exp(-coef)", 
                paste("lower .", round(100 * conf.int, 2), sep = ""), paste("upper .", 
                  round(100 * conf.int, 2), sep = "")))
        }
        rval$conf.int <- tmp
    }
    
    rval$varnames <- object$varnames
    rval$sig2 <- object$sig2
    rval$sig2var <- object$sig2var
    rval$rho <- object$rho
    rval$rhovar <- object$rhovar
    rval$nu <- object$nu
    rval$nuvar <- object$nuvar
    
    
    
    
    class(rval) <- "summary.parmfrailty"
    rval
}
