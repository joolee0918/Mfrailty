# ========================================================= 
# M-step for margins and gaussian copula 
# =========================================================


#----- Two stage  ------------------------------



Q1_ts <- function(p, n, Elog2u) {
    
    sig2 = exp(p)
    res = -n * log(sqrt(sig2)) - n * 1/8 * sig2 - (Elog2u)/(2 * sig2)
    return(-res)
    
}

Q1_g_ts <- function(p, n, Eu, Elogu) {
    
    sig2 = exp(p)
    res = -n * lgamma(1/sig2) - n * 1/sig2 * log(sig2) + (1/sig2 - 1) * (Elogu) - 
        Eu/sig2
    return(-res)
    
}


Q2_ts <- function(p, n, Eu, Elogu, Elog2u, Elogu1logu2, sig2, J, com) {
    
    res = 0
    K = J * (J - 1)/2
    for (i in 1:J) {
        tmp = 0
        tind <- which(com == i, arr.ind = T)
        for (j in 1:(J - 1)) tmp = tmp + Elog2u[[tind[j, 2]]][tind[j, 1]]
        
        res = res - n * log(sqrt(sig2[i])) - n * 1/8 * sig2[i] - tmp/(4 * sig2[i])
        
    }
    
    for (k in 1:K) {
        k1 = com[1, k]
        k2 = com[2, k]
        res = res + (-n/2 * 1/2 * log(1 - p[k]^2) - 1/2 * p[k]^2/(2 * (1 - p[k]^2)) * 
            ((Elog2u[[k]][1])/sig2[k1] + (Elogu[[k]][1]) + n * sig2[k1]/4 + (Elog2u[[k]][2])/sig2[k2] + 
                (Elogu[[k]][2]) + n * sig2[k2]/4) + 1/2 * p[k]/(1 - p[k]^2) * ((Elogu1logu2[[k]])/sqrt(sig2[k1] * 
            sig2[k2]) + sqrt(sig2[k2]) * (Elogu[[k]][1])/(2 * sqrt(sig2[k1])) + sqrt(sig2[k1]) * 
            (Elogu[[k]][2])/(2 * sqrt(sig2[k2])) + n * sqrt(sig2[k1] * sig2[k2])/4))
        
    }
    
    return(-res)
    
}


#----- gamma margins and Gaussian copula --------------

Q2_g_ts <- function(p, n, Lam, ww, vv, N, nev_id, parallel, ncore, dc, F1, f1, F2, 
    f2, J, com) {
    
    rho <- p
    
    K = J * (J - 1)/2
    
    res = 0
    
    newdc <- lognewdc <- list()
    
    if (any(rho <= -1) | any(rho >= 1)) {
        res <- NA
    } else {
        
        for (k in 1:K) {
            
            newdc[[k]] = (dCopula(cbind(eval(F1[[k]], list(x = vv[, 1])), eval(F2[[k]], 
                list(x = vv[, 2]))), normalCopula(rho[k])) * eval(f1[[k]], list(x = vv[, 
                1])) * eval(f2[[k]], list(x = vv[, 2])))
            
            lognewdc[[k]] <- ifelse(newdc[[k]] == 0, -double.xmax, log(newdc[[k]]))
            
            res = res + Ef_cal(f = lognewdc[[k]], cumsum1 = Lam[[com[1, k]]], cumsum2 = Lam[[com[2, 
                k]]], dc = dc[[k]], ww1 = ww[, 1], ww2 = ww[, 2], vv1 = vv[, 1], 
                vv2 = vv[, 2], nsample = n, N = N, n1 = nev_id[[com[1, k]]], n2 = nev_id[[com[2, 
                  k]]], parallel = parallel, ncore = ncore)
            
            
        }
        
    }
    return(-res)
    
}




#----- Pairwise M-step  -----------------------------------
#----- Lognormal margins and Gaussian copula --------------

Q2 <- function(p, n, Eu, Elogu, Elog2u, Elogu1logu2, J, com) {
    
    sig2 <- exp(p[1:J])
    rho <- p[-c(1:J)]
    
    K = J * (J - 1)/2
    
    
    res = 0
    
    for (i in 1:J) {
        tmp = 0
        tind <- which(com == i, arr.ind = T)
        for (j in 1:(J - 1)) tmp = tmp + Elog2u[[tind[j, 2]]][tind[j, 1]]
        
        res = res - n * log(sqrt(sig2[i])) - n * 1/8 * sig2[i] - tmp/(4 * sig2[i])
        
    }
    for (k in 1:K) {
        k1 = com[1, k]
        k2 = com[2, k]
        res = res + (-n/2 * 1/2 * log(1 - rho[k]^2) - 1/2 * rho[k]^2/(2 * (1 - rho[k]^2)) * 
            ((Elog2u[[k]][1])/sig2[k1] + (Elogu[[k]][1]) + n * sig2[k1]/4 + (Elog2u[[k]][2])/sig2[k2] + 
                (Elogu[[k]][2]) + n * sig2[k2]/4) + 1/2 * rho[k]/(1 - rho[k]^2) * 
            ((Elogu1logu2[[k]])/sqrt(sig2[k1] * sig2[k2]) + sqrt(sig2[k2]) * (Elogu[[k]][1])/(2 * 
                sqrt(sig2[k1])) + sqrt(sig2[k1]) * (Elogu[[k]][2])/(2 * sqrt(sig2[k2])) + 
                n * sqrt(sig2[k1] * sig2[k2])/4))
        
        
    }
    
    
    return(-res)
    
}



#----- general margins and Gaussian copula --------------

Q2_g <- function(p, n, Lam, ww, vv, N, nev_id, parallel, ncore, dc, J, com) {
    
    sig2 <- exp(p[1:J])
    rho <- p[-c(1:J)]
    
    K = J * (J - 1)/2
    
    res = 0
    
    newdc <- lognewdc <- list()
    
    if (any(rho <= -1) | any(rho >= 1)) {
        res <- NA
    } else {
        
        
        for (k in 1:K) {
            F1 <- copula:::asCall(paste0("p", "gamma"), list(shape = 1/sig2[com[1, 
                k]], scale = sig2[com[1, k]]))
            f1 <- copula:::asCall(paste0("d", "gamma"), list(shape = 1/sig2[com[1, 
                k]], scale = sig2[com[1, k]]))
            F2 <- copula:::asCall(paste0("p", "gamma"), list(shape = 1/sig2[com[2, 
                k]], scale = sig2[com[2, k]]))
            f2 <- copula:::asCall(paste0("d", "gamma"), list(shape = 1/sig2[com[2, 
                k]], scale = sig2[com[2, k]]))
            
            newdc[[k]] = (dCopula(cbind(eval(F1, list(x = vv[, 1])), eval(F2, list(x = vv[, 
                2]))), normalCopula(rho[k])) * eval(f1, list(x = vv[, 1])) * eval(f2, 
                list(x = vv[, 2])))
            
            lognewdc[[k]] <- ifelse(newdc[[k]] == 0, -double.xmax, log(newdc[[k]]))
            
            res = res + Ef_cal(f = lognewdc[[k]], cumsum1 = Lam[[com[1, k]]], cumsum2 = Lam[[com[2, 
                k]]], dc = dc[[k]], ww1 = ww[, 1], ww2 = ww[, 2], vv1 = vv[, 1], 
                vv2 = vv[, 2], nsample = n, N = N, n1 = nev_id[[com[1, k]]], n2 = nev_id[[com[2, 
                  k]]], parallel = parallel, ncore = ncore)
            
            
        }
    }
    
    return(-res)
    
}


