# ========================================================= 
# Propfile Likelihood with fixed beta, sigma, rho
# =========================================================


#----- Pairwise Likelihood  ------------------------------

parPL <- function(formula, data, frailty, theta, margins, vv, ww, nsample, N, nev_id, 
    num, idx, J, com, eventnames, parallel, control, dist, param) {
    
    
    K = J * (J - 1)/2
    
    mf <- Terms <- Y <- type <- X <- id <- offset <- offset.names <- list()
    
    
    for (i in 1:J) {
        mf[[i]] <- model.frame(formula[[i]], data[[i]])
        Terms[[i]] <- terms(mf[[i]])
    }
    
    
    ## Y, X, Frailty ID
    
    for (i in 1:J) {
        Y[[i]] <- model.extract(mf[[i]], "response")
        type[[i]] <- attr(Y[[i]], "type")
        if (type[[i]] != "right" && type[[i]] != "counting") 
            stop(paste("Cox model doesn't support \"", type[[i]], "\" survival data", 
                "for event", eventnames[i], sep = ""))
        
        X[[i]] <- model.matrix(formula[[i]], data[[i]])
        X[[i]] <- X[[i]][, -1, drop = FALSE]
        id[[i]] <- data[[i]][, frailty]
        offset[i] <- list(model.offset(mf[[i]]))
        offset.names[[i]] <- untangle.offset(Terms[[i]])
        if (is.null(offset[[i]]) | all(offset[[i]] == 0)) 
            offset[[i]] <- rep(0, nrow(mf[[i]])) else if (any(!is.finite(offset[[i]]))) 
            stop("offsets must be finite", "for event", eventnames[i], sep = "")
        if (sum(Y[[i]][, ncol(Y[[i]])]) == 0) {
            stop("No events", eventnames[i], sep = " ", "in the data")
        }
    }
    
    ncov <- sapply(1:J, function(i) ncol(X[[i]]))
    
    start <- 1
    beta <- alpha <- list()
    sig2 <- rep(0, J)
    
    for (i in 1:J) {
        beta[[i]] <- theta[start:(start + ncov[i] - 1)]
        alpha[[i]] <- exp(theta[(start + ncov[i]):(start + ncov[i] + 1)])
        sig2[i] <- theta[(start + ncov[i] + 2)]
        start <- start + ncov[i] + 3
    }
    rho <- theta[-c(1:(sum(ncov) + 2 * J + J))]
    
    
    shape <- scale <- hazards <- Hazards <- lam <- Lam <- list()
    
    for (i in 1:J) {
        shape[[i]] <- alpha[[i]][1]
        scale[[i]] <- alpha[[i]][2]
        
        if (dist[i] == "weibull") {
            hazards[[i]] <- eha::hweibull(Y[[i]][, 2], shape = shape[[i]], scale = scale[[i]])
            Hazards[[i]] <- eha::Hweibull(Y[[i]][, 2], shape = shape[[i]], scale = scale[[i]]) - 
                eha::Hweibull(Y[[i]][, 1], shape = shape[[i]], scale = scale[[i]])
        } else if (dist[i] == "loglogistic") {
            hazards[[i]] <- eha::hllogis(Y[[i]][, 2], shape = shape[[i]], scale = scale[[i]])
            Hazards[[i]] <- eha::Hllogis(Y[[i]][, 2], shape = shape[i], scale = scale[i]) - 
                eha::Hllogis(Y[[i]][, 1], shape = shape[i], scale = scale[i])
        } else if (dist[i] == "lognormal") {
            sdlog <- 1/shape[[i]]
            meanlog <- log(scale[[i]])
            
            hazards[[i]] <- eha::hlnorm(Y[[i]][, 2], meanlog = meanlog[i], sdlog = sdlog[i])
            Hazards[[i]] <- eha::Hlnorm(Y[[i]][, 2], meanlog = meanlog[i], sdlog = sdlog[i]) - 
                eha::Hlnorm(Y[[i]][, 1], meanlog = meanlog[i], sdlog = sdlog[i])
        } else if (dist[i] == "ev") {
            hazards[[i]] <- eha::hEV(Y[[i]][, 2], shape = shape[i], scale = scale[i])
            Hazards[[i]] <- eha::HEV(Y[[i]][, 2], shape = shape[i], scale = scale[i]) - 
                eha::HEV(Y[[i]][, 1], shape = shape[i], scale = scale[i])
        } else {
            stop(paste(dist, " is not availalble"))
        }
        lam[[i]] <- (hazards[[i]] * exp(X[[i]] %*% beta[[i]] + offset[[i]]))^(Y[[i]][, 
            3])
        lam[[i]] <- aggregate(lam[[i]], list(id[[i]]), prod, drop = T, simplify = T)[, 
            -1]
        Lam[[i]] <- drop(rowsum(Hazards[[i]] * exp(X[[i]] %*% beta[[i]] + offset[[i]]), 
            id[[i]], reorder = FALSE))
    }
    
    
    logLH <- dc <- list()
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
        
        # Expectation
        logLH[[k]] = ObslogL(haz1 = lam[[com[1, k]]], haz2 = lam[[com[2, k]]], cumsum1 = Lam[[com[1, 
            k]]], cumsum2 = Lam[[com[2, k]]], dc = dc[[k]], ww1 = ww[, 1], ww2 = ww[, 
            2], vv1 = vv[, 1], vv2 = vv[, 2], nsample = nsample, N = N, n1 = nev_id[[com[1, 
            k]]], n2 = nev_id[[com[2, k]]], parallel = parallel, ncore = control$ncore)
    }
    
    res <- Reduce(`+`, logLH)/(J - 1)
    
    return(res)
    
}

