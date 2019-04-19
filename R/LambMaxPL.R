# ========================================================= 
# Propfile Likelihood with fixed beta, sigma, rho
# =========================================================


#----- Pairwise Likelihood  ------------------------------

LambMaxPL <- function(formula, data, frailty, theta, cumhaz, Lam, margins, vv, ww, 
    nsample, N, nev_id, num, idx, J, com, eventnames, parallel, control) {
    
    
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
    p1_beta = sapply(1:J, function(i) ncol(X[[i]]))
    beta <- list()
    sig2 <- rep(0, J)
    p1_start = 1
    p1_end = -1
    
    for (l in 1:J) {
        p1_end = p1_end + p1_beta[l] + 1
        beta[[l]] <- theta[(p1_start):(p1_end)]
        sig2[l] <- theta[(p1_end + 1)]
        p1_start = p1_start + p1_beta[l] + 1
    }
    
    
    
    rho <- theta[(1 + sum(p1_beta) + J):(length(theta))]
    
    oldlam <- oldLam <- list()
    newformula <- newmodel <- newhazards <- newcumhaz <- newlam <- newLam <- list()
    Ef <- dc <- list()
    
    oldcumhaz = cumhaz
    oldLam = Lam
    
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
    
    
    
    error = 1
    iter = 0
    
    while (error > control$tol.P & iter < control$maxiter) {
        
        for (k in 1:K) {
            # Expectation
            Ef[[k]] = EU_cal(cumsum1 = oldLam[[com[1, k]]], cumsum2 = oldLam[[com[2, 
                k]]], dc = dc[[k]], ww1 = ww[, 1], ww2 = ww[, 2], vv1 = vv[, 1], 
                vv2 = vv[, 2], nsample = nsample, N = N, n1 = nev_id[[com[1, k]]], 
                n2 = nev_id[[com[2, k]]], parallel = parallel, ncore = control$ncore)
            
        }
        
        ## Updated Cumulative Hazard
        
        
        for (i in 1:J) {
            tind <- which(com == i, arr.ind = T)
            tmp = rep(0, nsample)
            tmpdata = data[[i]]
            for (j in 1:(J - 1)) tmp = tmp + Ef[[tind[j, 2]]][, tind[j, 1]]
            
            EU <- log(rep(tmp/(J - 1), num[[i]]))
            # print(EU)
            tmpdata$EU <- EU
            tmpdata$XEU <- X[[i]] %*% beta[[i]] + EU
            
            newformula[[i]] <- update.formula(formula[[i]], as.formula(~+offset(XEU)))
            if (!all(offset[[i]] == 0)) 
                newformula[[i]] <- update.formula(newformula[[i]], paste("~.+", paste(offset.names[[i]], 
                  collapse = "+")))
            newmodel[[i]] <- coxph(newformula[[i]], data = tmpdata, model = T, method = "breslow")
            
            scores <- exp(tmpdata$XEU + offset[[i]])
            stratum <- rep(1, NROW(Y[[i]]))
            newhazards[[i]] <- eha:::getHaz(Y[[i]], stratum, scores)[[1]]
            newhazards[[i]] <- rbind(c(0, 0), newhazards[[i]])
            newcumhaz[[i]] <- cumsum(newhazards[[i]][, 2])
            newlam[[i]] <- (newhazards[[i]][idx[[i]], 2] * exp(X[[i]] %*% beta[[i]] + 
                offset[[i]]))^{
                Y[[i]][, 3]
            }
            newlam[[i]] <- aggregate(newlam[[i]], list(id[[i]]), prod, drop = T, 
                simlify = T)[, -1]
            
            newLam[[i]] <- drop(rowsum(predict(newmodel[[i]], type = "expected")/exp(tmpdata$EU), 
                id[[i]], reorder = FALSE))
            
            
        }
        
        error = max(abs(unlist(lapply(1:J, function(j) newcumhaz[[j]] - oldcumhaz[[j]]))))
        
        iter = iter + 1
        # print(c(error, iter, sum(plogLH)))
        oldcumhaz = newcumhaz
    }
    
    logLH <- list()
    for (k in 1:K) {
        # Expectation
        logLH[[k]] = ObslogL(haz1 = newlam[[com[1, k]]], haz2 = newlam[[com[2, k]]], 
            cumsum1 = newLam[[com[1, k]]], cumsum2 = newLam[[com[2, k]]], dc = dc[[k]], 
            ww1 = ww[, 1], ww2 = ww[, 2], vv1 = vv[, 1], vv2 = vv[, 2], nsample = nsample, 
            N = N, n1 = nev_id[[com[1, k]]], n2 = nev_id[[com[2, k]]], parallel = parallel, 
            ncore = control$ncore)
    }
    
    plogLH <- Reduce(`+`, logLH)/(J - 1)
    
    
    res <- plogLH
    
    return(res)
    
}


#----- Two-stage Pairwise Likelihood  ---------------------------

LambMaxPL_ts <- function(formula, data, frailty, theta, cumhaz, Lam, margins, v, 
    w, nsample, N, nev_id, num, idx, eventnames, parallel, control) {
    
    
    mf <- model.frame(formula, data)
    Terms <- terms(mf)
    
    ## Y, X, frailty ID
    
    
    Y <- model.extract(mf, "response")
    type <- attr(Y, "type")
    if (type != "right" && type != "counting") 
        stop(paste("Cox model doesn't support \"", type, "\" survival data", "for event", 
            eventnames, sep = ""))
    
    X <- model.matrix(formula, data)
    X <- X[, -1, drop = FALSE]
    id <- data[, frailty]
    offset <- model.offset(mf)
    offset.names <- untangle.offset(Terms)
    if (is.null(offset) | all(offset == 0)) 
        offset <- rep(0, nrow(mf)) else if (any(!is.finite(offset))) 
        stop("offsets must be finite", "for event", eventnames, sep = "")
    if (sum(Y[, ncol(Y)]) == 0) 
        stop("No events", eventnames, sep = " ", "in the data")
    
    p1_beta = ncol(X)
    
    beta = theta[1:p1_beta]
    sig2 = theta[-c(1:p1_beta)]
    
    
    oldcumhaz = cumhaz
    oldLam = Lam
    
    error = 1
    iter = 0
    while (error > control$tol.P & iter < control$maxiter) {
        
        # Expectation of U
        
        Ef = EU_cal_ts(cumsum = oldLam, sig2 = sig2, margins = margins, w = w, v = v, 
            nsample = nsample, N = N, n = nev_id, parallel = parallel, ncore = control$ncore)
        
        ## Updated Cumulative Hazard
        
        tmpdata = data
        
        EU <- log(rep(Ef, num))
        # print(EU)
        tmpdata$EU <- EU
        tmpdata$XEU <- X %*% beta + EU
        
        newformula <- update.formula(formula, as.formula(~+offset(XEU)))
        if (!all(offset == 0)) 
            newformula <- update.formula(newformula, paste("~.+", paste(offset.names, 
                collapse = "+")))
        
        newmodel <- coxph(newformula, data = tmpdata, model = T, method = "breslow")
        
        scores <- exp(tmpdata$XEU + offset)
        stratum <- rep(1, NROW(Y))
        newhazards <- eha:::getHaz(Y, stratum, scores)[[1]]
        newhazards <- rbind(c(0, 0), newhazards)
        newcumhaz <- cumsum(newhazards[, 2])
        newlam <- (newhazards[idx, 2] * exp(X %*% beta + offset))^{
            Y[, 3]
        }
        newlam <- aggregate(newlam, list(id), prod, drop = T, simplify = T)[, -1]
        
        newLam <- drop(rowsum(predict(newmodel, type = "expected")/exp(tmpdata$EU), 
            id, reorder = FALSE))
        
        
        error = max(abs(newcumhaz - oldcumhaz))
        
        iter = iter + 1
        oldcumhaz = newcumhaz
    }
    
    # Expectation
    logLH = ObslogL_ts(haz = newlam, cumsum = newLam, sig2 = sig2, margins = margins, 
        w = w, v = v, nsample = nsample, N = N, n = nev_id, parallel = parallel, 
        ncore = control$ncore)
    
    
    res <- logLH
    
    return(res)
}

