# ==================================================================== 
# Estimation based on Two-stage Pairwise Likelihood with EM algorithm
# ====================================================================

tspairEM <- function(beta, cumhaz, sig2, rho, lam, Lam, formula, data, frailty, J, 
    paircop, margins, idx, num, control, w, v, ww, vv, nsample, N, com, nev_id, eventnames, 
    parallel) {
    
    mf <- Y <- type <- X <- id <- offset <- list()
    
    for (i in 1:J) {
        mf[[i]] <- model.frame(formula[[i]], data[[i]])
    }
    conv1 = rep(0, J)
    conv2 = 0
    
    
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
    Ef <- MLH <- PLH <- oldmlogLH <- newmlogLH <- list()
    
    conv1 = rep(0, J)
    iter1 = rep(0, J)
    
    
    K = J * (J - 1)/2
    oldbeta = beta
    oldcumhaz = cumhaz
    oldrho = rho
    oldsig2 = sig2
    oldlam = lam
    oldLam = Lam
    oldpsi = c(unlist(oldbeta), oldsig2, unlist(oldcumhaz))
    
    newsig2 = rep(0, J)
    
    ### First Stage: Estimate beta, cumhaz, sig2 ###
    
    for (i in 1:J) {
        oldmlogLH[[i]] = 0
        error1 = error2 = 1
        oldpsi = c(unlist(oldbeta[[i]]), oldsig2[i], unlist(oldcumhaz[[i]]))
        
        while (error1 > control$tol.P | iter1[i] < control$maxiter) {
            
            MLH[[i]] = ObslogL_ts(haz = oldlam[[i]], cumsum = oldLam[[i]], sig2 = oldsig2[i], 
                margins = margins[i], w = w, v = v, nsample = nsample, N = sqrt(N), 
                n = nev_id[[i]], parallel = parallel, ncore = control$ncore)
            
            if (anyNA(MLH[[i]])) {
                warnings("fail to calculate likelihood")
                conv1[i] = 1
                break
            }
            
            # Observed Liklihood
            
            newmlogLH[[i]] = sum(MLH[[i]])
            
            error2 = abs(newmlogLH[[i]] - oldmlogLH[[i]])/(abs(oldmlogLH[[i]]) + 
                .Machine$double.eps * 2)
            
            if (error2 < control$tol.L) 
                break
            
            #----- E-step -----------------------------------------------------------------#
            # Expectation of U
            Ef[[i]] = EU_cal_ts(cumsum = oldLam[[i]], sig2 = oldsig2[i], margins = margins[i], 
                w = w, v = v, nsample = nsample, N = sqrt(N), n = nev_id[[i]], parallel = parallel, 
                ncore = control$ncore)
            
            if (anyNA(Ef[[i]])) {
                warnings("fail to calculate E[U]")
                conv1[i] = 1
                break
            }
            
            #----- M-step -----------------------------------------------------------------#
            ## Updated beta, Cumulative Hazard
            
            tmpdata = data[[i]]
            
            tmpdata$EU <- log(rep(Ef[[i]], num[[i]]))
            newformula[[i]] <- update.formula(formula[[i]], as.formula(~. + offset(EU)))
            newmodel[[i]] <- coxph(newformula[[i]], data = tmpdata, model = T, method = "breslow")
            
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
            
            ## Updated sig2
            if (margins[i] == 1) {
                Elog2u_ts = Ef_cal_ts(f = (log(v))^2, cumsum = oldLam[[i]], sig2 = oldsig2[i], 
                  margins = margins[i], w = w, v = v, nsample = nsample, N = sqrt(N), 
                  n = nev_id[[i]], parallel = parallel, ncore = control$ncore)
                
                suppressWarnings(M.sig2 <- optimize(Q1_ts, c(-1e+08, 10), n = nsample, Elog2u = Elog2u_ts))
                
            } else {
                Elogu_ts = Ef_cal_ts(f = (log(v)), cumsum = oldLam[[i]], sig2 = oldsig2[i], 
                  margins = margins[i], w = w, v = v, nsample = nsample, N = sqrt(N), 
                  n = nev_id[[i]], parallel = parallel, ncore = control$ncore)
                
                M.sig2 <- optimize(Q1_g_ts, c(-1e+08, 10), n = nsample, Elogu = Elogu_ts, 
                  Eu = sum(Ef[[i]]))
                
            }
            
            newsig2[i] <- exp(M.sig2$minimum)
            
            newpsi = c(unlist(newbeta[[i]]), newsig2[i], unlist(newcumhaz[[i]]))
            error1 = max(abs(newpsi - oldpsi))
            iter1[i] = iter1[i] + 1
            oldbeta[[i]] = newbeta[[i]]
            oldlam[[i]] = newlam[[i]]
            oldcumhaz[[i]] = newcumhaz[[i]]
            oldLam[[i]] = newLam[[i]]
            oldsig2[i] = newsig2[i]
            oldpsi = newpsi
            oldmlogLH[[i]] = newmlogLH[[i]]
            
            
        }
        
        if (iter1[i] == control$maxiter) {
            warning(paste("Stage 1: did not converge in ", control$maxiter, " iterations."))
            conv1[i] = 1
            break
        }
    }
    
    
    if (any(conv1 == 1) | paircop == FALSE) {
        
        omega <- list()
        omega$beta <- newbeta
        omega$haz <- newhazards
        omega$cumhaz <- newcumhaz
        omega$lam <- newlam
        omega$Lam <- newLam
        omega$sig2 <- newsig2
        omega$loglik <- newmlogLH
        omega$conv <- c(conv1)
        omega$iter <- c(iter1)
        
    } else {
        ##### Second Stage #### Estimate rho
        
        oldplogLH = 0
        F1 <- f1 <- F2 <- f2 <- dc <- list()
        error3 = error4 = 1
        iter2 = 0
        for (k in 1:K) {
            if (margins[com[1, k]] == 1) {
                F1[[k]] <- copula:::asCall(paste0("p", "lnorm"), list(meanlog = -newsig2[com[1, 
                  k]]/2, sdlog = sqrt(newsig2[com[1, k]])))
                f1[[k]] <- copula:::asCall(paste0("d", "lnorm"), list(meanlog = -newsig2[com[1, 
                  k]]/2, sdlog = sqrt(newsig2[com[1, k]])))
            } else {
                F1[[k]] <- copula:::asCall(paste0("p", "gamma"), list(shape = 1/newsig2[com[1, 
                  k]], scale = newsig2[com[1, k]]))
                f1[[k]] <- copula:::asCall(paste0("d", "gamma"), list(shape = 1/newsig2[com[1, 
                  k]], scale = newsig2[com[1, k]]))
            }
            if (margins[com[2, k]] == 1) {
                F2[[k]] <- copula:::asCall(paste0("p", "lnorm"), list(meanlog = -newsig2[com[2, 
                  k]]/2, sdlog = sqrt(newsig2[com[2, k]])))
                f2[[k]] <- copula:::asCall(paste0("d", "lnorm"), list(meanlog = -newsig2[com[2, 
                  k]]/2, sdlog = sqrt(newsig2[com[2, k]])))
            } else {
                F2[[k]] <- copula:::asCall(paste0("p", "gamma"), list(shape = 1/newsig2[com[2, 
                  k]], scale = newsig2[com[2, k]]))
                f2[[k]] <- copula:::asCall(paste0("d", "gamma"), list(shape = 1/newsig2[com[2, 
                  k]], scale = newsig2[com[2, k]]))
            }
        }
        
        nlower = rep(lower, K)
        nupper = rep(upper, K)
        
        while (error3 > control$tol.P | iter2 < control$maxiter) {
            
            for (k in 1:K) {
                dc[[k]] = (dCopula(cbind(eval(F1[[k]], list(x = vv[, 1])), eval(F2[[k]], 
                  list(x = vv[, 2]))), normalCopula(oldrho[k])) * eval(f1[[k]], list(x = vv[, 
                  1])) * eval(f2[[k]], list(x = vv[, 2])))
                
            }
            
            if (all(margins == 1)) {
                Elogu <- Elog2u <- Elogu1logu2 <- Eu <- list()
                for (k in 1:K) {
                  
                  Eu[[k]] = colSums(EU_cal(cumsum1 = oldLam[[com[1, k]]], cumsum2 = oldLam[[com[2, 
                    k]]], dc = dc[[k]], ww1 = ww[, 1], ww2 = ww[, 2], vv1 = vv[, 
                    1], vv2 = vv[, 2], nsample = nsample, N = N, n1 = nev_id[[com[1, 
                    k]]], n2 = nev_id[[com[2, k]]], parallel = parallel, ncore = control$ncore))
                  
                  
                  Elogu[[k]] = Ef1f2_cal(f1 = log(vv[, 1]), f2 = log(vv[, 2]), cumsum1 = oldLam[[com[1, 
                    k]]], cumsum2 = oldLam[[com[2, k]]], dc = dc[[k]], ww1 = ww[, 
                    1], ww2 = ww[, 2], vv1 = vv[, 1], vv2 = vv[, 2], nsample = nsample, 
                    N = N, n1 = nev_id[[com[1, k]]], n2 = nev_id[[com[2, k]]], parallel = parallel, 
                    ncore = control$ncore)
                  
                  Elog2u[[k]] = Ef1f2_cal(f1 = (log(vv[, 1]))^2, f2 = (log(vv[, 2]))^2, 
                    cumsum1 = oldLam[[com[1, k]]], cumsum2 = oldLam[[com[2, k]]], 
                    dc = dc[[k]], ww1 = ww[, 1], ww2 = ww[, 2], vv1 = vv[, 1], vv2 = vv[, 
                      2], nsample = nsample, N = N, n1 = nev_id[[com[1, k]]], n2 = nev_id[[com[2, 
                      k]]], parallel = parallel, ncore = control$ncore)
                  
                  Elogu1logu2[[k]] = Ef_cal(f = log(vv[, 1]) * log(vv[, 2]), cumsum1 = oldLam[[com[1, 
                    k]]], cumsum2 = oldLam[[com[2, k]]], dc = dc[[k]], ww1 = ww[, 
                    1], ww2 = ww[, 2], vv1 = vv[, 1], vv2 = vv[, 2], nsample = nsample, 
                    N = N, n1 = nev_id[[com[1, k]]], n2 = nev_id[[com[2, k]]], parallel = parallel, 
                    ncore = control$ncore)
                  
                }
                
            }
            
            for (k in 1:K) {
                # Observed Liklihood
                PLH[[k]] = ObslogL(haz1 = newlam[[com[1, k]]], haz2 = newlam[[com[2, 
                  k]]], cumsum1 = newLam[[com[1, k]]], cumsum2 = newLam[[com[2, k]]], 
                  dc = dc[[k]], ww1 = ww[, 1], ww2 = ww[, 2], vv1 = vv[, 1], vv2 = vv[, 
                    2], nsample = nsample, N = N, n1 = nev_id[[com[1, k]]], n2 = nev_id[[com[2, 
                    k]]], parallel = parallel, ncore = control$ncore)
                
            }
            
            
            if (any(sapply(1:K, function(k) anyNA(PLH[[k]])))) {
                warning("fail to calculate likelihood")
                conv2 = 1
                break
            }
            
            newplogLH = sum(Reduce(`+`, PLH)/(J - 1))
            
            error4 = abs(newplogLH - oldplogLH)/(abs(oldplogLH) + .Machine$double.eps * 
                2)
            if (error4 < control$tol.L) 
                break
            
            if ( all(margins == 1)) {
                M.rho <- optim(par = oldrho, fn = Q2_ts, n = nsample, Eu = Eu, Elogu = Elogu, 
                  Elog2u = Elog2u, Elogu1logu2 = Elogu1logu2, sig2 = newsig2, J = J, 
                  com = com)
                newrho <- M.rho$par
            } else {
                M.phi <- optim(par = c(oldrho), fn = Q2_g_ts, n = nsample, Lam = oldLam, 
                  ww = ww, vv = vv, N = N, nev_id = nev_id, parallel = parallel, 
                  ncore = control$ncore, F1 = F1, f1 = f1, F2 = F2, f2 = f2, dc = dc, 
                  J = J, com = com)
                newrho <- M.phi$par
            }
            rm(M.rho)
            
            error3 = max(abs(newrho - oldrho))
            iter2 = iter2 + 1
            oldrho = newrho
            oldplogLH = newplogLH
            
            if (any((newrho - 2 * control$delta) <= -1) | any((newrho + 2 * control$delta) >= 
                1)) {
                message("Dependence parameter", " for event ", eventnames[((newrho - 
                  2 * control$delta) <= -1) | ((newrho + 2 * control$delta) >= 1)], 
                  " reaches the boundary", collapse = " ")
                conv2 <- 2
                break
            }
            
        }
        
        if (iter2 == control$maxiter) {
            warning(paste("Stage 2: did not converge in ", control$maxiter, " iterations."))
            conv2 <- 1
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
        omega$conv <- c(conv1, conv2)
        omega$iter <- c(iter1, iter2)
        
    }
    return(omega)
    
}
