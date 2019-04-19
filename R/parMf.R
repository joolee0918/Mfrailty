# ========================================================= 
# Maximize composite likelihood for parametric model
# =========================================================

#----- Pairwise M-step  -----------------------------------

mPL <- function(par, J, Y, X, offset, dist, margins, ww, vv, nsample, N, id, nev_id, 
    com, parallel, ncore) {
    
    K = J * (J - 1)/2
    ncov = sapply(1:J, function(i) ncol(X[[i]]))
    
    start <- 1
    beta <- alpha <- list()
    sig2 <- rep(0, J)
    
    for (i in 1:J) {
        beta[[i]] <- par[start:(start + ncov[i] - 1)]
        alpha[[i]] <- exp(par[(start + ncov[i]):(start + ncov[i] + 1)])
        sig2[i] <- exp(par[(start + ncov[i] + 2)])
        start <- start + ncov[i] + 3
    }
    rho <- par[-c(1:(sum(ncov) + 2 * J + J))]
    
    if (any(rho <= -1) | any(rho >= 1)) {
        res <- NA
    } else {
        newhazards <- newHazards <- newlam <- newLam <- list()
        # Expectation
        for (i in 1:J) {
            shape = alpha[[i]][1]
            scale = alpha[[i]][2]
            
            
            if (dist[i] == "weibull") {
                newhazards[[i]] <- eha::hweibull(Y[[i]][, 2], shape = shape, scale = scale)
                newHazards[[i]] <- eha::Hweibull(Y[[i]][, 2], shape = shape, scale = scale) - 
                  eha::Hweibull(Y[[i]][, 1], shape = shape, scale = scale)
            } else if (dist[i] == "loglogistic") {
                newhazards[[i]] <- eha::hllogis(Y[[i]][, 2], shape = shape, scale = scale)
                newHazards[[i]] <- eha::Hllogis(Y[[i]][, 2], shape = shape, scale = scale) - 
                  eha::Hllogis(Y[[i]][, 1], shape = shape, scale = scale)
            } else if (dist[i] == "lognormal") {
                sdlog <- 1/shape
                meanlog <- log(scale)
                
                newhazards[[i]] <- eha::hlnorm(Y[[i]][, 2], meanlog = meanlog, sdlog = sdlog)
                newHazards[[i]] <- eha::Hlnorm(Y[[i]][, 2], meanlog = meanlog, sdlog = sdlog) - 
                  eha::Hlnorm(Y[[i]][, 1], meanlog = meanlog, sdlog = sdlog)
            } else if (dist[i] == "ev") {
                newhazards[[i]] <- eha::hEV(Y[[i]][, 2], shape = shape, scale = scale)
                newHazards[[i]] <- eha::HEV(Y[[i]][, 2], shape = shape, scale = scale) - 
                  eha::HEV(Y[[i]][, 1], shape = shape, scale = scale)
                
            } else {
                stop(paste(dist, " is not availalble"))
            }
            
            
            newlam[[i]] <- (newhazards[[i]] * exp(X[[i]] %*% beta[[i]] + offset[[i]]))^(Y[[i]][, 
                3])
            newlam[[i]] <- aggregate(newlam[[i]], list(id[[i]]), prod, drop = T, 
                simplify = T)[, -1]
            newLam[[i]] <- drop(rowsum(newHazards[[i]] * exp(X[[i]] %*% beta[[i]] + 
                offset[[i]]), id[[i]], reorder = FALSE))
        }
        
        PLH <- dc <- list()
        
        for (k in 1:K) {
            if (margins[com[1, k]] == 1) {
                F1 <- copula:::asCall(paste0("p", "lnorm"), list(meanlog = -sig2[com[1, 
                  k]]/2, sdlog = sqrt(sig2[com[1, k]])))
                f1 <- copula:::asCall(paste0("d", "lnorm"), list(meanlog = -sig2[com[1, 
                  k]]/2, sdlog = sqrt(sig2[com[1, k]])))
            } else {
                F1 <- copula:::asCall(paste0("p", "gamma"), list(shape = 1/sig2[com[1, 
                  k]], scale = sig2[com[1, k]]))
                f1 <- copula:::asCall(paste0("d", "gamma"), list(shape = 1/sig2[com[1, 
                  k]], scale = sig2[com[1, k]]))
            }
            if (margins[com[2, k]] == 1) {
                F2 <- copula:::asCall(paste0("p", "lnorm"), list(meanlog = -sig2[com[2, 
                  k]]/2, sdlog = sqrt(sig2[com[2, k]])))
                f2 <- copula:::asCall(paste0("d", "lnorm"), list(meanlog = -sig2[com[2, 
                  k]]/2, sdlog = sqrt(sig2[com[2, k]])))
            } else {
                F2 <- copula:::asCall(paste0("p", "gamma"), list(shape = 1/sig2[com[2, 
                  k]], scale = sig2[com[2, k]]))
                f2 <- copula:::asCall(paste0("d", "gamma"), list(shape = 1/sig2[com[2, 
                  k]], scale = sig2[com[2, k]]))
            }
            
            dc[[k]] = (dCopula(cbind(eval(F1, list(x = vv[, 1])), eval(F2, list(x = vv[, 
                2]))), normalCopula(rho[k])) * eval(f1, list(x = vv[, 1])) * eval(f2, 
                list(x = vv[, 2])))
        }
        
        
        for (k in 1:K) {
            # Observed Liklihood
            PLH[[k]] = ObslogL(haz1 = newlam[[com[1, k]]], haz2 = newlam[[com[2, 
                k]]], cumsum1 = newLam[[com[1, k]]], cumsum2 = newLam[[com[2, k]]], 
                dc = dc[[k]], ww1 = ww[, 1], ww2 = ww[, 2], vv1 = vv[, 1], vv2 = vv[, 
                  2], nsample = nsample, N = N, n1 = nev_id[[com[1, k]]], n2 = nev_id[[com[2, 
                  k]]], parallel = parallel, ncore = ncore)
        }
        
        res = -sum(Reduce(`+`, PLH)/(J - 1))
    }
    return(res)
    
}

