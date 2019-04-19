# ========================================================= 
# Dependence modeling for a multivariate mixed-Poisson model via copulas
# =========================================================



#' Fitting a semiparametric multivariate mixed-Poisson model
#'
#' Fitting a semiparametric multivariate mixed-Poisson model via a Gaussian copula with a log-normal marginal distribution of random effects
#'
#' @param formula a formula object with an obect of the type \code{Surv} on the left side and the terms on the right side. Note that \code{+strata()} is not supported.
#' @param data a 'data.frame' which has variables in 'formula'
#' @param frailty a charhacter string specifying a group
#' @param paircop a logical value: if \code{TRUE}, a dependent pairwise likelihood is used; otherwise, an independent mixed-Poissom model is used. 
#' @param ntype the number of types of events.
#' @param twostage a logical value: if \code{TRUE}, two-stage procedure is used. 
#' @param model logical value: if \code{TRUE}, the model frame is returned.
#' @param model.matrix logical value: if \code{TRUE}, the x matirx is returned.
#' @param control an object of class specifying control options created by \code{\link{controlList}}.
#' @param init a list specifying inital values of a vector of variances of random effects and dependence parameters. Default initial value for a variance of random effect is 1 and 0 for the dependence parameter.
#' @param eventnames a vector of character string speficying event names
#' @param parallel logical value: if \code{TRUE}, it supports parallel execution.
#' @param ncore the number of cores if \code{parallel = TRUE}.

#' @return an object of class \code{emmfrailty} representing the fit.
#' \item{beta.coefficients}{a list of each type of events of the estimated regressoin coefficients.}
#' \item{sig2}{a vector of the estimated variances of random effects.}
#' \item{rho}{a vector of the estimated dependence parameters from a Gaussian copula function.}
#' \item{ktau}{a vector of the estimated Kendall's tau from a Gaussian copula function.}
#' \item{basecmhaz}{a list of each type of events of the estiamted cumulative baseline hazards.}
#' \item{logLik}{a value of log likelihood at the final values of the parameters.}
#' \item{iter}{the number of iterations used.}
#' \item{conv}{an integer code for the convergence. 0 indicates successful convergence, 1 indicates a failure to convergence, and 2 indicates the estimated dependence parameter reaches the boundary. }
#' \item{var}{a variance-covariance matrix of the estimates.}
#' \item{betavar}{a list of each type of the variances of the regression coefficient estimates.}
#' \item{sig2var}{a vector of the variances of random effects variances estimates.}
#' \item{rhovar}{a vector of the dependence parameter estimates.}
#' \item{varnames}{a list of event names and variable names.}
#' \item{control}{a list of control arguments used.}
#' \item{n}{the number of sample size.}
#' \item{nevent}{a vector of the number of each types of events.}
#' \item{neventtype}{the number of types of events.}
#' \item{nknot}{the number nodes used in Gaussian-quadrature.}
#' \item{margins}{a vector of character strings specyting marginal distributions of random effects.}
#' 
#' The object will also contain the following: two_stage, copula, model, call, optionally x, and model.
#


#' @importFrom stats model.frame model.extract model.matrix model.offset aggregate update.formula optim optimize printCoefmat rexp qlnorm rnorm pchisq pnorm predict plnorm dlnorm pgamma dgamma
#' @import survival 
#' @importFrom copula normalCopula rCopula iTau tau dMvdc mvdc optimMeth getTheta dCopula
#' @importFrom utils combn
#' @importFrom Matrix bdiag
#' @importFrom parallel detectCores

#' @export
emmfrailty <- function(formula, data, frailty = NULL, paircop = TRUE, margins = "lnorm", ntype, 
    twostage = FALSE, model = FALSE, model.matrix = FALSE, control = list(), init = list(), 
    eventnames = NULL, parallel = FALSE, ncore = NULL, ...) {
    
    
    if (is.null(ntype)) 
        stop("number of evnet types is missing")
    J = ntype
    K = ntype * (ntype - 1)/2
    com <- combn(1:J, 2)
    
    control <- controlList(control)
    
    if (is.null(eventnames)) 
        eventnames = seq(1:J)
    if (parallel == TRUE) {
        if (is.null(ncore)) 
            control$ncore = parallel::detectCores(logical = TRUE)/2 else control$ncore = ncore
    } else control$ncore = 1
    
    
    Call <- match.call()
    
    if (missing(formula) | missing(data)) 
        stop("Missing arguments")
    
    if(length(margins) != J) margins <- rep(margins, J)
    cmargins <- margins
    margins <- ifelse(margins == "lnorm", 1, 2)
    
    
    
    mf <- Y <- type <- X <- id <- offset <- offset.names <- list()
    
    for (i in 1:J) {
        mf[[i]] <- model.frame(formula[[i]], data[[i]])
     }
    
    ## Y, X, Frailty ID, offset
    
    for (i in 1:J) {
        Y[[i]] <- model.extract(mf[[i]], "response")
        type[[i]] <- attr(Y[[i]], "type")
        if (type[[i]] != "right" && type[[i]] != "counting") 
            stop(paste("Cox model doesn't support \"", type[[i]], "\" survival data", 
                "for event", eventnames[i], sep = ""))
        
        X[[i]] <- model.matrix(formula[[i]], data[[i]], contrast = NULL)
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
    
    ### number of events in ID
    order_id <- nev_id <- num <- list()
    nsample <- length(unique(id[[1]]))
    
    for (i in 1:J) {
        order_id[[i]] <- match(id[[i]], unique(id[[i]]))
        nev_id[[i]] <- as.numeric(rowsum(Y[[i]][, 3], order_id[[i]], reorder = FALSE))  # nevent per frailty
        num[[i]] <- as.numeric(table(id[[i]]))
    }
    
    # Gaussian Legendre node & weight --------------------------------
    
    gauss.quad <- gauleg.f(control$nknot, 0, 1)
    w <- gauss.quad$weight
    v <- -log(gauss.quad$location)
    
    ww <- expand.grid(w, w)
    vv <- expand.grid(v, v)
    N = control$nknot^2
    
    
    ############### EM algorithm ###################################################
    
    #----- Initial parameters ----------------------------------------#
    
    mcox <- newmodel <- newlog1 <- newformula <- list()
    oldcumhaz <- oldbeta <- oldlam <- oldLam <- oldhazards <- list()
    idx <- list()
    
    for (i in 1:J) {
        mcox[[i]] <- coxph(formula[[i]], data[[i]], model = T, method = "breslow")
        oldbeta[[i]] <- mcox[[i]]$coefficients
        scores <- exp(X[[i]] %*% oldbeta[[i]])
        stratum <- rep(1, NROW(Y[[i]]))
        oldhazards[[i]] <- eha:::getHaz(Y[[i]], stratum, scores)[[1]]
        oldhazards[[i]] <- rbind(c(0, 0), oldhazards[[i]])
        oldcumhaz[[i]] <- cumsum(oldhazards[[i]][, 2])
        
    }
    
    
    if (is.null(init$sig2)) {
        oldsig2 <- rep(1, J)
    } else oldsig2 <- init$sig2
    
    if (is.null(init$rho)) {
        oldrho <- rep(0, K)
    } else oldrho <- init$rho
    
    # Index for event time -----------------------
    for (i in 1:J) {
        idx[[i]] = findInterval(Y[[i]][, 2], vec = oldhazards[[i]][, 1])
    }
    
    # Hazard and cumulative hazard ---------------
    for (i in 1:J) {
        oldlam[[i]] <- (oldhazards[[i]][idx[[i]], 2] * as.vector(exp(X[[i]] %*% oldbeta[[i]])))^{
            Y[[i]][, 3]
        }
        oldlam[[i]] <- aggregate(oldlam[[i]], list(id[[i]]), prod, simplify = TRUE)[, 
            -1]
        oldLam[[i]] <- predict(mcox[[i]], type = "expected", collapse = id[[i]])
    }
    
    
    # EM Algorithm --------------------------------
    
    if (twostage == FALSE) {
        if (paircop == TRUE) {
            fit <- pairEM(beta = oldbeta, cumhaz = oldcumhaz, sig2 = oldsig2, 
                rho = oldrho, lam = oldlam, Lam = oldLam, formula = formula, 
                data = data, frailty = frailty, J = J, margins = margins, idx = idx, 
                num = num, control = control, ww = ww, vv = vv, nsample = nsample, 
                N = N, com = com, nev_id = nev_id, eventnames = eventnames, parallel = parallel)
        } else {
            fit <- tspairEM(beta = oldbeta, cumhaz = oldcumhaz, sig2 = oldsig2, rho = oldrho, 
                lam = oldlam, Lam = oldLam, formula = formula, data = data, frailty = frailty, 
                J = J, paircop = paircop, margins = margins, idx = idx, num = num, 
                control = control, w = w, v = v, ww = ww, vv = vv, nsample = nsample, 
                N = N, com = com, nev_id = nev_id, eventnames = eventnames, parallel = parallel)
        }
        
    } else {
        fit <- tspairEM(beta = oldbeta, cumhaz = oldcumhaz, sig2 = oldsig2, rho = oldrho, 
            lam = oldlam, Lam = oldLam, formula = formula, data = data, frailty = frailty, 
            J = J, paircop = paircop, margins = margins, idx = idx, num = num, control = control, 
            w = w, v = v, ww = ww, vv = vv, nsample = nsample, N = N, com = com, 
            nev_id = nev_id, eventnames = eventnames, parallel = parallel)
    }
    
    if (sum(fit$conv) > 0) {
        result <- list(fail = fit)
        message("Fail to fit")
    } else {
        
        # Estimates -------------------------------
        
        est_beta = fit$beta
        names(est_beta) = eventnames
        est_cumhaz = fit$cumhaz
        names(est_cumhaz) = eventnames
        est_haz = fit$haz
        names(est_haz) = eventnames
        est_sig2 = fit$sig2
        est_rho = fit$rho
        # est_ktau = fit$ktau
        est_lam = fit$lam
        names(est_lam) = eventnames
        est_Lam = fit$Lam
        names(est_Lam) = eventnames
        
        
        # Variance ---------------------------------- Use profile likelihood method
        # ------------------------------------------
        if (twostage == FALSE) {
            if (paircop == TRUE) {
                V <- Var(beta = est_beta, sig2 = est_sig2, rho = est_rho, 
                  lam = est_lam, Lam = est_Lam, cumhaz = est_cumhaz, formula = formula, 
                  data = data, frailty = frailty, J = J, margins = margins, idx = idx, 
                  num = num, control = control, ww = ww, vv = vv, nsample = nsample, 
                  N = N, com = com, nev_id = nev_id, eventnames = eventnames, parallel = parallel)
            } else {
                V <- tsVar(beta = fit$beta, sig2 = fit$sig2, rho = fit$rho,  
                  lam = fit$lam, Lam = fit$Lam, cumhaz = fit$cumhaz, formula = formula, 
                  data = data, frailty = frailty, J = J, paircop = paircop, margins = margins, 
                  idx = idx, num = num, control = control, w = w, v = v, ww = ww, 
                  vv = vv, nsample = nsample, N = N, com = com, nev_id = nev_id, 
                  eventnames = eventnames, parallel = parallel)
            }
        } else {
            V <- tsVar(beta = est_beta, sig2 = est_sig2, rho = est_rho,  
                lam = est_lam, Lam = est_Lam, cumhaz = est_cumhaz, formula = formula, 
                data = data, frailty = frailty, J = J, paircop = paircop, margins = margins, 
                idx = idx, num = num, control = control, w = w, v = v, ww = ww, vv = vv, 
                nsample = nsample, N = N, com = com, nev_id = nev_id, eventnames = eventnames, 
                parallel = parallel)
        }
        
        
        p1_beta = unlist(lapply(1:J, function(i) length(est_beta[[i]])))
        p1 = sum(p1_beta)
        p2 = length(est_sig2)
        if (paircop == TRUE) 
            p3 = length(c(est_rho))
        
        if (anyNA(V)) {
            result <- list(fail = fit)
            message("no standard error is availble")
        } else {
            betaV <- list()
            sig2V <- rep(0, J)
            p1_start = 1
            p1_end = -1
            diagV = diag(V)
            
            for (i in 1:J) {
                p1_end = p1_end + p1_beta[i] + 1
                betaV[[i]] <- diagV[(p1_start):(p1_end)]
                sig2V[i] <- diagV[(p1_end + 1)]
                p1_start = p1_start + p1_beta[i] + 1
            }
            sig2V <- sig2V
            if (paircop == TRUE) 
                rhoV <- diagV[(1 + p1 + p2):(p1 + p2 + p3)]
            
            newktau <- rep(0, K)
            
            if (paircop == TRUE) {
                for (i in 1:K) {
                  tmpcopula <- asCall0("normalCopula", est_rho[i])
                  newktau[i] <- tau(eval(tmpcopula))
                }
            }
            
            varNames <- list()
            varNames$eventnames <- eventnames
            varNames$beta.names <- lapply(1:J, function(i) paste("event", eventnames[i], 
                names(est_beta[[i]]), sep = "."))
            varNames$sig2.names <- paste("frailty", margins, "sd", sep = ".")
            varNames$rho.names <- paircop
            
            result <- list()
            result$omega <- c(unlist(est_beta), est_sig2, est_rho)
            result$beta.coefficients <- est_beta
            result$sig2 <- est_sig2
            if (paircop == TRUE) 
                result$rho <- est_rho
            result$ktau <- newktau
            result$basecmhaz <- lapply(1:J, function(i) cbind(est_cumhaz[[i]], est_haz[[i]][, 
                1]))
            result$logLik <- fit$loglik
            result$iter <- fit$iter
            result$var <- V
            result$betavar <- betaV
            result$sig2var <- sig2V
            if (paircop == TRUE) 
                result$rhovar <- rhoV
            result$varnames <- varNames
            result$convergence <- fit$conv
            result$control <- control
            result$n <- nsample
            result$nevent <- unlist(lapply(1:J, function(i) sum(nev_id[[i]])))
            result$neventtype <- J
            result$nknot <- control$nkont
            result$margins <- cmargins
            if (paircop == TRUE) 
                result$copula <- "normalCopula"
            result$two_stage <- twostage
            
            if (model.matrix) 
                result$x <- X
            if (model) 
                result$model <- mf
            result$call <- Call
            
            
        }
    }
    attr(result, "class") <- "emmfrailty"
    
    return(result)
    
}
