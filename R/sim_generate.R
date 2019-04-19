# ======================================================= 
# Generate 3 types of recurrent data Random Effects: Multivariate lognormal distribution with Normal
# copula Event Times: Truncated weibull distribution
# =======================================================



TWeiRandom.f <- function(tt, uu, lam, alp, tau) {
    v <- runif(1, min = 0, max = 1)
    
    nomin <- (lam * (tt^alp) * uu) - log(1 - v)
    denom <- lam * uu
    term <- (nomin/denom)^(1/alp)
    term <- ifelse(is.na(term), tau, term)
    return(term)
}

#' @export
generatedata.f <- function(nsample, J, copula, margins, tau, probC, sig21, lam1, 
    alp1, beta11, beta12, sig22, lam2, alp2, beta21, beta22, sig23, lam3, alp3, beta31, 
    beta32, rho, meanw, dispstr = NULL) {
    
    getdata.f <- function(id, trt, W, tau, lam, beta1, beta2, alp, uu) {
        
        lam <- lam * exp(beta1 * trt + beta2 * W)
        
        cur.t <- TWeiRandom.f(tt = 0, uu = uu, lam = lam, alp = alp, tau = tau)
        if (cur.t >= tau) {
            estart <- 0
            estop <- tau
            estatus <- 0
            enum <- 1
            nlen <- 0
        } else {
            vec.t <- cur.t
            while (cur.t < tau) {
                pre.t <- cur.t
                cur.t <- TWeiRandom.f(tt = pre.t, uu = uu, lam = lam, alp = alp, 
                  tau = tau)
                vec.t <- c(vec.t, cur.t)
            }
            vec.t <- vec.t[vec.t < tau]
            nlen <- length(vec.t)
            
            id <- rep(id, nlen + 1)
            estart <- c(0, vec.t)
            estop <- c(vec.t, tau)
            estatus <- c(rep(1, nlen), 0)
            enum <- 1:(nlen + 1)
        }
        
        tmp <- data.frame(id = id, estart = estart, estop = estop, estatus = estatus, 
            enum = enum, n = nlen, tau = tau, trt = trt, W = W)
        return(tmp)
    }
    
    if (probC == 0) {
        CC <- rep(tau, nsample)
    } else {
        CC <- rexp(nsample, rate = ((-1) * log(1 - probC)))
        CC <- ifelse(CC > tau, tau, CC)
    }
    
    
    Copula = asCall2(copula, param = rho, dim = J, dispstr = dispstr)
    
    paramMargins <- list()
    
    if (margins[1] == "lnorm") 
        paramMargins[[1]] <- list(meanlog = -sig21/2, sdlog = sqrt(sig21)) else paramMargins[[1]] <- list(shape = 1/sig21, scale = sig21)
    if (margins[2] == "lnorm") 
        paramMargins[[2]] <- list(meanlog = -sig22/2, sdlog = sqrt(sig22)) else paramMargins[[2]] <- list(shape = 1/sig22, scale = sig22)
    if (margins[3] == "lnorm") 
        paramMargins[[3]] <- list(meanlog = -sig23/2, sdlog = sqrt(sig23)) else paramMargins[[3]] <- list(shape = 1/sig23, scale = sig23)
    
    simcop <- rCopula(nsample, eval(Copula))
    U1 <- eval(copula:::asCall(paste0("q", margins[1]), paramMargins[[1]]), list(x = simcop[, 
        1]))
    U2 <- eval(copula:::asCall(paste0("q", margins[2]), paramMargins[[2]]), list(x = simcop[, 
        2]))
    U3 <- eval(copula:::asCall(paste0("q", margins[3]), paramMargins[[3]]), list(x = simcop[, 
        3]))
    
    # print(c(sum(log(U1)^2), sum(log(U2)^2), sum(log(U3)^2)))
    # print(c(sum(log(U1)*log(U2)), sum(log(U1)*log(U3)),sum(log(U2)*log(U3))))
    
    trt <- c(rep(0, nsample/2), rep(1, nsample/2))
    W <- rnorm(nsample, meanw, 1)
    
    event1 <- lapply(1:nsample, function(i) getdata.f(id = i, trt = trt[i], W = W[i], 
        tau = CC[i], lam = lam1, beta1 = beta11, beta2 = beta12, alp = alp1, uu = U1[i]))
    event2 <- lapply(1:nsample, function(i) getdata.f(id = i, trt = trt[i], W = W[i], 
        tau = CC[i], lam = lam2, beta1 = beta21, beta2 = beta22, alp = alp2, uu = U2[i]))
    event3 <- lapply(1:nsample, function(i) getdata.f(id = i, trt = trt[i], W = W[i], 
        tau = CC[i], lam = lam3, beta1 = beta31, beta2 = beta32, alp = alp3, uu = U3[i]))
    
    event1 <- do.call("rbind", event1)
    event2 <- do.call("rbind", event2)
    event3 <- do.call("rbind", event3)
    
    return(list(event1, event2, event3))
}
