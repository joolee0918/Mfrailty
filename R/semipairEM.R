# ========================================================= 
# Estimation based on Pairwise Likelihood with EM algorithm
# =========================================================

pairEM <- function(beta, cumhaz, sig2, rho, lam, Lam, formula, data, frailty, J, 
    margins, idx, num, control, ww, vv, nsample, N, com, nev_id, eventnames, parallel) {
    
    mf <- Y <- type <- X <- id <- offset <- list()
    
    for (i in 1:J) {
        mf[[i]] <- model.frame(formula[[i]], data[[i]])
    }
    
    K = J * (J - 1)/2
    conv = 0
    
    ## Y, X, frailty ID
    
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
        if (is.null(offset[[i]]) | all(offset[[i]] == 0)) 
            offset[[i]] <- rep(0, nrow(mf[[i]])) else if (any(!is.finite(offset[[i]]))) 
            stop("offsets must be finite", "for event", eventnames[i], sep = "")
        if (sum(Y[[i]][, ncol(Y[[i]])]) == 0) {
            stop("No events", eventnames[i], sep = " ", "in the data")
        }
    }
    
    
    tmpcopula <- asCall0("normalCopula", NULL)
    cop.param <- getTheta(eval(tmpcopula), freeOnly = TRUE, attr = TRUE)
    lower <- attr(cop.param, "param.lowbnd")
    upper <- attr(cop.param, "param.upbnd")
        
    #----- E-step -----------------------------------------------------------------#
    newformula <- newmodel <- newlam <- newcumhaz <- newhazards <- newbeta <- newLam <- list()
    Ef <- PLH <- dc <- list()
    
    K = J * (J - 1)/2
    oldbeta = beta
    oldcumhaz = cumhaz
    oldsig2 = sig2
    oldrho = rho
    
    oldlam = lam
    oldLam = Lam
    
    oldplogLH = 0
    error1 = error2 = 1
    iter = 0
    
    oldomega = c(unlist(oldbeta), oldsig2, oldrho, unlist(oldcumhaz))
    
    while (error1 > control$tol.P | iter < control$maxiter) {
        
        for (k in 1:K) {
            if (margins[com[1, k]] == 1) {
                F1 <- copula:::asCall(paste0("p", "lnorm"), list(meanlog = -oldsig2[com[1, 
                  k]]/2, sdlog = sqrt(oldsig2[com[1, k]])))
                f1 <- copula:::asCall(paste0("d", "lnorm"), list(meanlog = -oldsig2[com[1, 
                  k]]/2, sdlog = sqrt(oldsig2[com[1, k]])))
            } else {
                F1 <- copula:::asCall(paste0("p", "gamma"), list(shape = 1/oldsig2[com[1, 
                  k]], scale = oldsig2[com[1, k]]))
                f1 <- copula:::asCall(paste0("d", "gamma"), list(shape = 1/oldsig2[com[1, 
                  k]], scale = oldsig2[com[1, k]]))
            }
            if (margins[com[2, k]] == 1) {
                F2 <- copula:::asCall(paste0("p", "lnorm"), list(meanlog = -oldsig2[com[2, 
                  k]]/2, sdlog = sqrt(oldsig2[com[2, k]])))
                f2 <- copula:::asCall(paste0("d", "lnorm"), list(meanlog = -oldsig2[com[2, 
                  k]]/2, sdlog = sqrt(oldsig2[com[2, k]])))
            } else {
                F2 <- copula:::asCall(paste0("p", "gamma"), list(shape = 1/oldsig2[com[2, 
                  k]], scale = oldsig2[com[2, k]]))
                f2 <- copula:::asCall(paste0("d", "gamma"), list(shape = 1/oldsig2[com[2, 
                  k]], scale = oldsig2[com[2, k]]))
            }
            
            dc[[k]] = (dCopula(cbind(eval(F1, list(x = vv[, 1])), eval(F2, list(x = vv[, 
                2]))), normalCopula(oldrho[k])) * eval(f1, list(x = vv[, 1])) * eval(f2, 
                list(x = vv[, 2])))
            
            # Observed Liklihood
            PLH[[k]] = ObslogL(haz1 = oldlam[[com[1, k]]], haz2 = oldlam[[com[2, 
                k]]], cumsum1 = oldLam[[com[1, k]]], cumsum2 = oldLam[[com[2, k]]], 
                dc = dc[[k]], ww1 = ww[, 1], ww2 = ww[, 2], vv1 = vv[, 1], vv2 = vv[, 
                  2], nsample = nsample, N = N, n1 = nev_id[[com[1, k]]], n2 = nev_id[[com[2, 
                  k]]], parallel = parallel, ncore = control$ncore)
        }
        
        if (any(sapply(1:K, function(k) anyNA(PLH[[k]])))) {
            warnings("fail to calculate likelihood")
            conv2 = 1
            break
        }
        
        newplogLH = sum(Reduce(`+`, PLH)/(J - 1))
        error2 = abs(newplogLH - oldplogLH)/(abs(oldplogLH) + .Machine$double.eps * 
            2)
        
        # if (error1 < control$tol.P | error2 < control$tol.L)
        if (error2 < control$tol.L) 
            break
        
        #----- E-step -----------------------------------------------------------------#
        # Expectation of U
        
        for (k in 1:K) {
            
            Ef[[k]] = EU_cal(cumsum1 = oldLam[[com[1, k]]], cumsum2 = oldLam[[com[2, 
                k]]], dc = dc[[k]], ww1 = ww[, 1], ww2 = ww[, 2], vv1 = vv[, 1], 
                vv2 = vv[, 2], nsample = nsample, N = N, n1 = nev_id[[com[1, k]]], 
                n2 = nev_id[[com[2, k]]], parallel = parallel, ncore = control$ncore)
            
        }
        
        if ( all(margins == 1)) {
            
            Elogu <- Elog2u <- Elogu1logu2 <- Eu <- list()
            
            for (k in 1:K) {
                
                Elogu[[k]] = Ef1f2_cal(f1 = log(vv[, 1]), f2 = log(vv[, 2]), cumsum1 = oldLam[[com[1, 
                  k]]], cumsum2 = oldLam[[com[2, k]]], dc = dc[[k]], ww1 = ww[, 1], 
                  ww2 = ww[, 2], vv1 = vv[, 1], vv2 = vv[, 2], nsample = nsample, 
                  N = N, n1 = nev_id[[com[1, k]]], n2 = nev_id[[com[2, k]]], parallel = parallel, 
                  ncore = control$ncore)
                
                Elog2u[[k]] = Ef1f2_cal(f1 = (log(vv[, 1]))^2, f2 = (log(vv[, 2]))^2, 
                  cumsum1 = oldLam[[com[1, k]]], cumsum2 = oldLam[[com[2, k]]], dc = dc[[k]], 
                  ww1 = ww[, 1], ww2 = ww[, 2], vv1 = vv[, 1], vv2 = vv[, 2], nsample = nsample, 
                  N = N, n1 = nev_id[[com[1, k]]], n2 = nev_id[[com[2, k]]], parallel = parallel, 
                  ncore = control$ncore)
                
                Elogu1logu2[[k]] = Ef_cal(f = log(vv[, 1]) * log(vv[, 2]), cumsum1 = oldLam[[com[1, 
                  k]]], cumsum2 = oldLam[[com[2, k]]], dc = dc[[k]], ww1 = ww[, 1], 
                  ww2 = ww[, 2], vv1 = vv[, 1], vv2 = vv[, 2], nsample = nsample, 
                  N = N, n1 = nev_id[[com[1, k]]], n2 = nev_id[[com[2, k]]], parallel = parallel, 
                  ncore = control$ncore)
                
                Eu[[k]] = colSums(Ef[[k]])
                
            }
            
        }
        
        
        if (any(sapply(1:K, function(k) anyNA(Ef[[k]])))) {
            warnings("fail to calculate E[U]")
            conv2 = 1
            break
        }
        
        #----- M-step -----------------------------------------------------------------#
        ## Updated beta and Cumulative Hazard
        
        for (i in 1:J) {
            tind <- which(com == i, arr.ind = T)
            tmp = rep(0, nsample)
            for (j in 1:(J - 1)) tmp = tmp + Ef[[tind[j, 2]]][, tind[j, 1]]
            
            tmpdata <- data[[i]]
            
            tmpdata$EU <- log(rep(tmp/(J - 1), num[[i]]))
            newformula[[i]] <- stats::update.formula(formula[[i]], as.formula(~. + 
                offset(EU)))
            # message(print(head(tmpdata)))
            newmodel[[i]] <- survival::coxph(newformula[[i]], data = tmpdata, model = T, 
                method = "breslow")
            
            newbeta[[i]] <- newmodel[[i]]$coefficients
            scores <- exp(tmpdata$EU + X[[i]] %*% newbeta[[i]] + offset[[i]])
            stratum <- rep(1, NROW(Y[[i]]))
            newhazards[[i]] <- eha:::getHaz(Y[[i]], stratum, scores)[[1]]
            newhazards[[i]] <- rbind(c(0, 0), newhazards[[i]])
            newcumhaz[[i]] <- cumsum(newhazards[[i]][, 2])
            newlam[[i]] <- (newhazards[[i]][idx[[i]], 2] * exp(X[[i]] %*% newbeta[[i]] + 
                offset[[i]]))^(Y[[i]][, 3])
            newlam[[i]] <- aggregate(newlam[[i]], list(id[[i]]), prod, drop = T, 
                simplify = T)[, -1]
            newLam[[i]] <- drop(rowsum(predict(newmodel[[i]], type = "expected")/exp(tmpdata$EU), 
                id[[i]], reorder = FALSE))
            
        }
        
        ## Updated sigma, rho
        
        if (all(margins == 1)) {
            M.phi <- optim(par = c(log(oldsig2), oldrho), fn = Q2, n = nsample, Eu = Eu, 
                Elogu = Elogu, Elog2u = Elog2u, Elogu1logu2 = Elogu1logu2, J = J, 
                com = com)
            newsig2 <- exp(M.phi$par[1:J])
            newrho <- M.phi$par[-c(1:J)]
        } else {
            M.phi <- optim(par = c(log(oldsig2), oldrho), fn = Q2_g, n = nsample, 
                Lam = oldLam, ww = ww, vv = vv, N = N, nev_id = nev_id, parallel = parallel, 
                ncore = control$ncore, dc = dc, J = J, com = com)
            newsig2 <- exp(M.phi$par[1:J])
            newrho <- M.phi$par[-c(1:J)]
        }
        
        rm(M.phi)
        
        
        newomega <- c(unlist(newbeta), newsig2, newrho, unlist(newcumhaz))
        
        error1 <- max(abs(newomega - oldomega)/(abs(oldomega) + .Machine$double.eps * 
            2))
        iter <- iter + 1
        oldbeta = newbeta
        oldcumhaz = newcumhaz
        oldlam = newlam
        oldLam = newLam
        oldsig2 = newsig2
        oldrho = newrho
        oldomega = newomega
        oldplogLH = newplogLH
        
        
        
        if (any((newrho - 2 * control$delta) <= -1) | any((newrho + 2 * control$delta) >= 
            1)) {
            message("Dependence parameter", " for event ", eventnames[((newrho - 
                2 * control$delta) <= -1) | ((newrho + 2 * control$delta) >= 1)], 
                " reaches the boundary", collapse = " ")
            conv <- 2
            break
        }
        
        
    }
    
    if (iter == control$maxiter) {
        warnings(paste("did not converge in ", control$maxiter, " iterations."))
        conv <- 1
    }
    omega <- list()
    
    omega$beta <- newbeta
    omega$haz <- newhazards
    omega$cumhaz <- newcumhaz
    omega$lam <- newlam
    omega$Lam <- newLam
    omega$sig2 <- newsig2
    omega$rho <- newrho
    omega$loglik <- newplogLH
    omega$conv <- conv
    omega$iter <- iter
    
    
    return(omega)
    
}
