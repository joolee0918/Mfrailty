# ========================================================= Variance and
# Covariance Calculation based on Profile LH
# =========================================================

#----- Pairwise Likelihood  -----------------------------------

Var <- function(beta, sig2, rho, lam, Lam, cumhaz, formula, data, frailty, 
    J, margins, idx, num, control, ww, vv, nsample, N, com, nev_id, eventnames, parallel) {
    
    K = J * (J - 1)/2
    
    
    theta1 = unlist(lapply(1:J, function(i) c(beta[[i]], sig2[i])))
    theta = c(theta1, rho)
    len = length(theta)
    
    delta = control$delta
    
    
    tmpLH <- LambMaxPL(formula = formula, data = data, frailty = frailty, theta = theta, 
        cumhaz = cumhaz, Lam = Lam, margins = margins, vv = vv, ww = ww, nsample = nsample, 
        N = N, nev_id = nev_id, idx = idx, num = num, J = J, com = com, eventnames = eventnames, 
        parallel = parallel, control = control)
    
    pl <- sum(tmpLH)/nsample
    pl_i <- tmpLH
    
    
    PLs1 <- matrix(0, len, len)
    for (i in 1:len) {
        for (j in i:len) {
            the1 <- theta
            the1[i] <- theta[i] + delta
            the1[j] <- the1[j] + delta
            
            tmpLH <- LambMaxPL(formula = formula, data = data, frailty = frailty, 
                theta = the1, cumhaz = cumhaz, Lam = Lam, margins = margins, vv = vv, 
                ww = ww, nsample = nsample, N = N, nev_id = nev_id, idx = idx, num = num, 
                J = J, com = com, eventnames = eventnames, parallel = parallel, control = control)
            
            PLs1[i, j] <- sum(tmpLH)/nsample
        }
    }
    
    
    PLs2 <- matrix(0, len, len)
    for (i in 1:len) {
        for (j in i:len) {
            the1 <- theta
            the1[i] <- theta[i] - delta
            the1[j] <- the1[j] + delta
            
            tmpLH <- LambMaxPL(formula = formula, data = data, frailty = frailty, 
                theta = the1, cumhaz = cumhaz, Lam = Lam, margins = margins, vv = vv, 
                ww = ww, nsample = nsample, N = N, nev_id = nev_id, idx = idx, num = num, 
                J = J, com = com, eventnames = eventnames, parallel = parallel, control = control)
            
            PLs2[i, j] <- sum(tmpLH)/nsample
        }
    }
    
    PLs3 <- matrix(0, len, len)
    for (i in 1:len) {
        for (j in i:len) {
            the1 <- theta
            the1[i] <- theta[i] + delta
            the1[j] <- the1[j] - delta
            
            tmpLH <- LambMaxPL(formula = formula, data = data, frailty = frailty, 
                theta = the1, cumhaz = cumhaz, Lam = Lam, margins = margins, vv = vv, 
                ww = ww, nsample = nsample, N = N, nev_id = nev_id, idx = idx, num = num, 
                J = J, com = com, eventnames = eventnames, parallel = parallel, control = control)
            
            PLs3[i, j] <- sum(tmpLH)/nsample
        }
    }
    
    PLs4 <- matrix(0, len, len)
    for (i in 1:len) {
        for (j in i:len) {
            the1 <- theta
            the1[i] <- theta[i] - delta
            the1[j] <- the1[j] - delta
            
            tmpLH <- LambMaxPL(formula = formula, data = data, frailty = frailty, 
                theta = the1, cumhaz = cumhaz, Lam = Lam, margins = margins, vv = vv, 
                ww = ww, nsample = nsample, N = N, nev_id = nev_id, idx = idx, num = num, 
                J = J, com = com, eventnames = eventnames, parallel = parallel, control = control)
            
            PLs4[i, j] <- sum(tmpLH)/nsample
        }
    }
    
    
    
    pl1 <- rep(0, len)
    pl1_i <- matrix(0, nrow = nsample, ncol = len)
    for (i in 1:len) {
        the1 <- theta
        the1[i] <- theta[i] + delta
        
        
        tmpLH <- LambMaxPL(formula = formula, data = data, frailty = frailty, theta = the1, 
            cumhaz = cumhaz, Lam = Lam, margins = margins, vv = vv, ww = ww, nsample = nsample, 
            N = N, nev_id = nev_id, idx = idx, num = num, J = J, com = com, eventnames = eventnames, 
            parallel = parallel, control = control)
        pl1[i] <- sum(tmpLH)/nsample
        pl1_i[, i] <- tmpLH
    }
    
    pl2 <- rep(0, len)
    pl2_i <- matrix(0, nrow = nsample, ncol = len)
    for (i in 1:len) {
        the1 <- theta
        the1[i] <- theta[i] - delta
        
        
        tmpLH <- LambMaxPL(formula = formula, data = data, frailty = frailty, theta = the1, 
            cumhaz = cumhaz, Lam = Lam, margins = margins, vv = vv, ww = ww, nsample = nsample, 
            N = N, nev_id = nev_id, idx = idx, num = num, J = J, com = com, eventnames = eventnames, 
            parallel = parallel, control = control)
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

#----- Two-stage Pairwise Likelihood  -------------------------------

tsVar <- function(beta, sig2, rho, lam, Lam, cumhaz, formula, data, frailty, 
    J, paircop, margins, idx, num, control, w, v, ww, vv, nsample, N, com, nev_id, 
    eventnames, parallel) {
    
    theta1 = lapply(1:J, function(i) c(beta[[i]], sig2[i]))
    theta2 = rho
    theta = c(unlist(theta1), theta2)
    
    len1 = lapply(1:J, function(i) length(theta1[[i]]))
    len2 = length(theta2)
    len = length(theta)
    
    p1_beta = lapply(1:J, function(i) length(beta[[i]]))
    p1 = sum(unlist(p1_beta))
    p2 = length(sig2)
    p3 = length(rho)
    
    delta = control$delta
    
    ml <- ml_i <- MLs1 <- MLs2 <- MLs3 <- MLs4 <- ml1 <- ml1_i <- ml2 <- ml2_i <- list()
    for (i in 1:J) {
        ### Observed likelihood
        
        tmpLH <- LambMaxPL_ts(formula = formula[[i]], data = data[[i]], frailty = frailty, 
            theta = theta1[[i]], cumhaz = cumhaz[[i]], Lam = Lam[[i]], margins = margins[i], 
            v = v, w = w, nsample = nsample, N = sqrt(N), nev_id = nev_id[[i]], idx = idx[[i]], 
            num = num[[i]], eventnames = eventnames[[i]], parallel = parallel, control = control)
        
        
        ml[[i]] <- sum(tmpLH)/nsample
        ml_i[[i]] <- tmpLH
        
        
        MLs1[[i]] <- matrix(0, len1[[i]], len1[[i]])
        for (j in 1:len1[[i]]) {
            for (k in j:len1[[i]]) {
                the1 <- theta1[[i]]
                the1[j] <- theta1[[i]][j] + delta
                the1[k] <- the1[k] + delta
                
                tmpLH <- LambMaxPL_ts(formula = formula[[i]], data = data[[i]], frailty = frailty, 
                  theta = the1, cumhaz = cumhaz[[i]], Lam = Lam[[i]], margins = margins[i], 
                  v = v, w = w, nsample = nsample, N = sqrt(N), nev_id = nev_id[[i]], 
                  idx = idx[[i]], num = num[[i]], eventnames = eventnames[[i]], parallel = parallel, 
                  control = control)
                
                MLs1[[i]][j, k] <- sum(tmpLH)/nsample
            }
        }
        
        MLs2[[i]] <- matrix(0, len1[[i]], len1[[i]])
        for (j in 1:len1[[i]]) {
            for (k in j:len1[[i]]) {
                the1 <- theta1[[i]]
                the1[j] <- theta1[[i]][j] - delta
                the1[k] <- the1[k] + delta
                
                tmpLH <- LambMaxPL_ts(formula = formula[[i]], data = data[[i]], frailty = frailty, 
                  theta = the1, cumhaz = cumhaz[[i]], Lam = Lam[[i]], margins = margins[i], 
                  v = v, w = w, nsample = nsample, N = sqrt(N), nev_id = nev_id[[i]], 
                  idx = idx[[i]], num = num[[i]], eventnames = eventnames[[i]], parallel = parallel, 
                  control = control)
                
                MLs2[[i]][j, k] <- sum(tmpLH)/nsample
            }
        }
        
        
        MLs3[[i]] <- matrix(0, len1[[i]], len1[[i]])
        for (j in 1:len1[[i]]) {
            for (k in j:len1[[i]]) {
                the1 <- theta1[[i]]
                the1[j] <- theta1[[i]][j] + delta
                the1[k] <- the1[k] - delta
                
                tmpLH <- LambMaxPL_ts(formula = formula[[i]], data = data[[i]], frailty = frailty, 
                  theta = the1, cumhaz = cumhaz[[i]], Lam = Lam[[i]], margins = margins[i], 
                  v = v, w = w, nsample = nsample, N = sqrt(N), nev_id = nev_id[[i]], 
                  idx = idx[[i]], num = num[[i]], eventnames = eventnames[[i]], parallel = parallel, 
                  control = control)
                
                MLs3[[i]][j, k] <- sum(tmpLH)/nsample
            }
        }
        
        
        MLs4[[i]] <- matrix(0, len1[[i]], len1[[i]])
        for (j in 1:len1[[i]]) {
            for (k in j:len1[[i]]) {
                the1 <- theta1[[i]]
                the1[j] <- theta1[[i]][j] - delta
                the1[k] <- the1[k] - delta
                
                tmpLH <- LambMaxPL_ts(formula = formula[[i]], data = data[[i]], frailty = frailty, 
                  theta = the1, cumhaz = cumhaz[[i]], Lam = Lam[[i]], margins = margins[i], 
                  v = v, w = w, nsample = nsample, N = sqrt(N), nev_id = nev_id[[i]], 
                  idx = idx[[i]], num = num[[i]], eventnames = eventnames[[i]], parallel = parallel, 
                  control = control)
                
                MLs4[[i]][j, k] <- sum(tmpLH)/nsample
            }
        }
        
        
        ml1[[i]] <- rep(0, len1[[i]])
        ml1_i[[i]] <- matrix(0, nrow = nsample, ncol = len1[[i]])
        for (j in 1:len1[[i]]) {
            the1 <- theta1[[i]]
            the1[j] <- theta1[[i]][j] + delta
            tmpLH <- LambMaxPL_ts(formula = formula[[i]], data = data[[i]], frailty = frailty, 
                theta = the1, cumhaz = cumhaz[[i]], Lam = Lam[[i]], margins = margins[i], 
                v = v, w = w, nsample = nsample, N = sqrt(N), nev_id = nev_id[[i]], 
                idx = idx[[i]], num = num[[i]], eventnames = eventnames[[i]], parallel = parallel, 
                control = control)
            ml1[[i]][j] <- sum(tmpLH)/nsample
            ml1_i[[i]][, j] <- tmpLH
        }
        
        ml2[[i]] <- rep(0, len1[[i]])
        ml2_i[[i]] <- matrix(0, nrow = nsample, ncol = len1[[i]])
        for (j in 1:len1[[i]]) {
            the1 <- theta1[[i]]
            the1[j] <- theta1[[i]][j] - delta
            tmpLH <- LambMaxPL_ts(formula = formula[[i]], data = data[[i]], frailty = frailty, 
                theta = the1, cumhaz = cumhaz[[i]], Lam = Lam[[i]], margins = margins[i], 
                v = v, w = w, nsample = nsample, N = sqrt(N), nev_id = nev_id[[i]], 
                idx = idx[[i]], num = num[[i]], eventnames = eventnames[[i]], parallel = parallel, 
                control = control)
            ml2[[i]][j] <- sum(tmpLH)/nsample
            ml2_i[[i]][, j] <- tmpLH
        }
        
    }
    
    
    if (paircop == TRUE) {
        K = J * (J - 1)/2
        
        theta1 = unlist(lapply(1:J, function(i) c(beta[[i]], sig2[i])))
        len11 = length(theta1)
        
        p1_beta = unlist(lapply(1:J, function(i) length(beta[[i]])))
        
        tmpLH <- LambMaxPL(formula = formula, data = data, frailty = frailty, theta = theta, 
            cumhaz = cumhaz, Lam = Lam, margins = margins, vv = vv, ww = ww, nsample = nsample, 
            N = N, nev_id = nev_id, idx = idx, num = num, J = J, com = com, eventnames = eventnames, 
            parallel = parallel, control = control)
        
        pl <- sum(tmpLH)/nsample
        pl_i <- tmpLH
        
        PLs1 <- matrix(0, len2, len)
        for (i in 1:len2) {
            for (j in 1:len) {
                the1 <- theta
                the1[len11 + i] <- theta[len11 + i] + delta
                the1[j] <- the1[j] + delta
                
                tmpLH <- LambMaxPL(formula = formula, data = data, frailty = frailty, 
                  theta = the1, cumhaz = cumhaz, Lam = Lam, margins = margins, vv = vv, 
                  ww = ww, nsample = nsample, N = N, nev_id = nev_id, idx = idx, 
                  num = num, J = J, com = com, eventnames = eventnames, parallel = parallel, 
                  control = control)
                PLs1[i, j] <- sum(tmpLH)/nsample
            }
        }
        PLs2 <- matrix(0, len2, len)
        for (i in 1:len2) {
            for (j in 1:len) {
                the1 <- theta
                the1[len11 + i] <- theta[len11 + i] - delta
                the1[j] <- the1[j] + delta
                
                tmpLH <- LambMaxPL(formula = formula, data = data, frailty = frailty, 
                  theta = the1, cumhaz = cumhaz, Lam = Lam, margins = margins, vv = vv, 
                  ww = ww, nsample = nsample, N = N, nev_id = nev_id, idx = idx, 
                  num = num, J = J, com = com, eventnames = eventnames, parallel = parallel, 
                  control = control)
                PLs2[i, j] <- sum(tmpLH)/nsample
            }
        }
        PLs3 <- matrix(0, len2, len)
        for (i in 1:len2) {
            for (j in 1:len) {
                the1 <- theta
                the1[len11 + i] <- theta[len11 + i] + delta
                the1[j] <- the1[j] - delta
                
                tmpLH <- LambMaxPL(formula = formula, data = data, frailty = frailty, 
                  theta = the1, cumhaz = cumhaz, Lam = Lam, margins = margins, vv = vv, 
                  ww = ww, nsample = nsample, N = N, nev_id = nev_id, idx = idx, 
                  num = num, J = J, com = com, eventnames = eventnames, parallel = parallel, 
                  control = control)
                PLs3[i, j] <- sum(tmpLH)/nsample
            }
        }
        
        PLs4 <- matrix(0, len2, len)
        for (i in 1:len2) {
            for (j in 1:len) {
                the1 <- theta
                the1[len11 + i] <- theta[len11 + i] - delta
                the1[j] <- the1[j] - delta
                
                tmpLH <- LambMaxPL(formula = formula, data = data, frailty = frailty, 
                  theta = the1, cumhaz = cumhaz, Lam = Lam, margins = margins, vv = vv, 
                  ww = ww, nsample = nsample, N = N, nev_id = nev_id, idx = idx, 
                  num = num, J = J, com = com, eventnames = eventnames, parallel = parallel, 
                  control = control)
                PLs4[i, j] <- sum(tmpLH)/nsample
            }
        }
        
        pl1 <- rep(0, len)
        pl1_i <- matrix(0, nrow = nsample, ncol = len)
        for (i in 1:len) {
            the1 <- theta
            the1[i] <- theta[i] + delta
            
            tmpLH <- LambMaxPL(formula = formula, data = data, frailty = frailty, 
                theta = the1, cumhaz = cumhaz, Lam = Lam, margins = margins, vv = vv, 
                ww = ww, nsample = nsample, N = N, nev_id = nev_id, idx = idx, num = num, 
                J = J, com = com, eventnames = eventnames, parallel = parallel, control = control)
            pl1[i] <- sum(tmpLH)/nsample
            pl1_i[, i] <- tmpLH
        }
        
        pl2 <- rep(0, len)
        pl2_i <- matrix(0, nrow = nsample, ncol = len)
        for (i in 1:len) {
            the1 <- theta
            the1[i] <- theta[i] - delta
            
            tmpLH <- LambMaxPL(formula = formula, data = data, frailty = frailty, 
                theta = the1, cumhaz = cumhaz, Lam = Lam, margins = margins, vv = vv, 
                ww = ww, nsample = nsample, N = N, nev_id = nev_id, idx = idx, num = num, 
                J = J, com = com, eventnames = eventnames, parallel = parallel, control = control)
            pl2[i] <- sum(tmpLH)/nsample
            pl2_i[, i] <- tmpLH
        }
        
    }
    
    
    I <- list()
    for (i in 1:J) {
        I[[i]] <- matrix(0, len1[[i]], len1[[i]])
        for (j in 1:len1[[i]]) {
            for (k in j:len1[[i]]) {
                if (j == k) 
                  I[[i]][j, j] <- -(ml1[[i]][j] - 2 * ml[[i]] + ml2[[i]][j])/(delta^2) else {
                  I[[i]][j, k] <- -(MLs1[[i]][j, k] - MLs2[[i]][j, k] - MLs3[[i]][j, 
                    k] + MLs4[[i]][j, k])/(4 * delta^2)
                  
                }
            }
        }
        I[[i]] <- I[[i]] + t(I[[i]]) - diag(diag(I[[i]]))
    }
    
    I <- as.matrix(bdiag(I))
    
    if (paircop == TRUE) {
        A2 <- matrix(0, len2, len)
        for (i in 1:len2) {
            for (j in 1:len) {
                if ((i + len11) == j) 
                  A2[i, j] <- -(pl1[j] - 2 * pl + pl2[j])/(delta^2) else {
                  A2[i, j] <- -(PLs1[i, j] - PLs2[i, j] - PLs3[i, j] + PLs4[i, j])/(4 * 
                    delta^2)
                }
            }
        }
        
        I <- rbind(cbind(I, matrix(0, nrow = len - len2, ncol = len2)), A2)
    }
    
    A <- try(solve(I), silent = TRUE)
    if (class(A) == "try-error") {
        warnings("Fail to inverse of information matrix")
        V <- NA
    } else {
        
        if (paircop == TRUE) {
            S1 <- list()
            
            for (i in 1:J) {
                S1[[i]] <- matrix(0, nrow = nsample, ncol = len1[[i]])
                for (j in 1:len1[[i]]) {
                  S1[[i]][, j] <- (ml1_i[[i]][, j] - ml2_i[[i]][, j])/(2 * delta)
                }
            }
            
            S2 <- matrix(0, nrow = nsample, ncol = len2)
            for (j in 1:len2) {
                S2[, j] <- (pl1_i[, len11 + j] - pl2_i[, len11 + j])/(2 * delta)
            }
            
            S = cbind(do.call(cbind, S1), S2)
            
            B <- t(S) %*% S/nsample
            
            
            V <- A %*% B %*% t(A)/nsample
        } else {
            V <- A/nsample
        }
    }
    return(V)
}
