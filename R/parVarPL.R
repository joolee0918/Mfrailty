# ========================================================= 
# Variance and Covariance Calculation based on Profile LH for a parametric model
# =========================================================

#----- Pairwise Likelihood  -----------------------------------

parVar <- function(beta, alpha, sig2, rho, nu = NULL, lam, Lam, formula, data, frailty, 
    J, margins, idx, num, control, ww, vv, nsample, N, com, nev_id, eventnames, parallel, 
    dist, param) {
    
    K = J * (J - 1)/2
    
    
    theta1 = unlist(lapply(1:J, function(i) c(beta[[i]], log(alpha[[i]]), sig2[i])))
    theta = c(theta1, rho)
    
    len = length(theta)
    
    
    delta = control$delta
    
    
    tmpLH <- parPL(formula = formula, data = data, frailty = frailty, theta = theta, 
        margins = margins, vv = vv, ww = ww, nsample = nsample, N = N, nev_id = nev_id, 
        idx = idx, num = num, J = J, com = com, eventnames = eventnames, parallel = parallel, 
        control = control, dist = dist, param = param)
    
    
    pl <- sum(tmpLH)/nsample
    pl_i <- tmpLH
    
    
    PLs1 <- matrix(0, len, len)
    for (i in 1:len) {
        for (j in i:len) {
            the1 <- theta
            the1[i] <- theta[i] + delta
            the1[j] <- the1[j] + delta
            
            tmpLH <- parPL(formula = formula, data = data, frailty = frailty, theta = the1, 
                margins = margins, vv = vv, ww = ww, nsample = nsample, N = N, nev_id = nev_id, 
                idx = idx, num = num, J = J, com = com, eventnames = eventnames, 
                parallel = parallel, control = control, dist = dist, param = param)
            
            
            PLs1[i, j] <- sum(tmpLH)/nsample
        }
    }
    
    
    PLs2 <- matrix(0, len, len)
    for (i in 1:len) {
        for (j in i:len) {
            the1 <- theta
            the1[i] <- theta[i] - delta
            the1[j] <- the1[j] + delta
            
            tmpLH <- parPL(formula = formula, data = data, frailty = frailty, theta = the1, 
                margins = margins, vv = vv, ww = ww, nsample = nsample, N = N, nev_id = nev_id, 
                idx = idx, num = num, J = J, com = com, eventnames = eventnames, 
                parallel = parallel, control = control, dist = dist, param = param)
            
            PLs2[i, j] <- sum(tmpLH)/nsample
        }
    }
    
    PLs3 <- matrix(0, len, len)
    for (i in 1:len) {
        for (j in i:len) {
            the1 <- theta
            the1[i] <- theta[i] + delta
            the1[j] <- the1[j] - delta
            
            tmpLH <- parPL(formula = formula, data = data, frailty = frailty, theta = the1, 
                margins = margins, vv = vv, ww = ww, nsample = nsample, N = N, nev_id = nev_id, 
                idx = idx, num = num, J = J, com = com, eventnames = eventnames, 
                parallel = parallel, control = control, dist = dist, param = param)
            
            PLs3[i, j] <- sum(tmpLH)/nsample
        }
    }
    
    PLs4 <- matrix(0, len, len)
    for (i in 1:len) {
        for (j in i:len) {
            the1 <- theta
            the1[i] <- theta[i] - delta
            the1[j] <- the1[j] - delta
            
            tmpLH <- parPL(formula = formula, data = data, frailty = frailty, theta = the1, 
                margins = margins, vv = vv, ww = ww, nsample = nsample, N = N, nev_id = nev_id, 
                idx = idx, num = num, J = J, com = com, eventnames = eventnames, 
                parallel = parallel, control = control, dist = dist, param = param)
            
            PLs4[i, j] <- sum(tmpLH)/nsample
        }
    }
    
    
    
    pl1 <- rep(0, len)
    pl1_i <- matrix(0, nrow = nsample, ncol = len)
    for (i in 1:len) {
        the1 <- theta
        the1[i] <- theta[i] + delta
        
        
        tmpLH <- parPL(formula = formula, data = data, frailty = frailty, theta = the1, 
            margins = margins, vv = vv, ww = ww, nsample = nsample, N = N, nev_id = nev_id, 
            idx = idx, num = num, J = J, com = com, eventnames = eventnames, parallel = parallel, 
            control = control, dist = dist, param = param)
        
        pl1[i] <- sum(tmpLH)/nsample
        pl1_i[, i] <- tmpLH
    }
    
    pl2 <- rep(0, len)
    pl2_i <- matrix(0, nrow = nsample, ncol = len)
    for (i in 1:len) {
        the1 <- theta
        the1[i] <- theta[i] - delta
        
        
        tmpLH <- parPL(formula = formula, data = data, frailty = frailty, theta = the1, 
            margins = margins, vv = vv, ww = ww, nsample = nsample, N = N, nev_id = nev_id, 
            idx = idx, num = num, J = J, com = com, eventnames = eventnames, parallel = parallel, 
            control = control, dist = dist, param = param)
        pl2[i] <- sum(tmpLH)/nsample
        pl2_i[, i] <- tmpLH
    }
    
    
    I <- matrix(0, len, len)
    for (i in 1:len) {
        for (j in i:len) {
            if (i == j) 
                I[i, i] <- -(pl1[i] - 2 * pl + pl2[i])/(delta^2) else {
                I[i, j] <- -(PLs1[i, j] - PLs2[i, j] - PLs3[i, j] + PLs4[i, j])/(4 * 
                  delta^2)
            }
        }
    }
    I <- I + t(I) - diag(diag(I))
    
    A <- try(solve(I), silent = TRUE)
    if (class(A) == "try-error") {
        warnings("Fail to inverse of information matrix")
        V <- NA
    } else {
        
        S <- matrix(0, nrow = nsample, ncol = len)
        for (j in 1:len) {
            S[, j] <- (pl1_i[, j] - pl2_i[, j])/(2 * delta)
        }
        
        
        B <- t(S) %*% S/nsample
        
        V <- A %*% B %*% t(A)/nsample
    }
    
    return(V)
}
