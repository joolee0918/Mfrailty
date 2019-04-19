# ========================================================= 
# Evaluate copula function Call function for marinal fraity distribution Transformation of Ktau
# to stablize estimation
# =========================================================


pairmvdc <- function(paircopula, margins, paramMargins, arginsIdentical = FALSE, 
    check = TRUE, fixupNames = TRUE) {
    K <- length(paircop)
    com <- combn(K, 2)
    
    pmvdc <- list()
    
    for (k in 1:K) {
        pmvdc[[k]] <- mvdc(eval(paircopula[[k]]), c(margins[com[1, k]], margins[com[2, 
            k]]), list(paramMargins[[com[1, k]]], paramMargins[[com[2, k]]]))
    }
    
    return(pmvdc)
}


asCall0 <- function(fun, param) {
    cc <- if (length(param) == 0) 
        quote(FUN()) else {
        as.call(c(quote(FUN), param))
    }
    
    cc[[1]] <- as.name(fun)
    cc
}

asCall2 <- function(fun, param, dim, dispstr) {
    cc <- if (length(param) == 0) 
        as.call(c(quote(FUN()), dim = dim, dispstr = dispstr)) else {
        as.call(c(quote(FUN), param = list(param), dim = dim, dispstr = dispstr))
    }
    
    cc[[1]] <- as.name(fun)
    cc
}

tranktau <- function(x) log((1 + x)/(1 - x))

itranktau <- function(x) (exp(x) - 1)/(exp(x) + 1)

ditranktau <- function(x) 2 * exp(x)/(exp(x) + 1)^2

