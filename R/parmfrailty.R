# ========================================================= 
# Parametric Dependence modeling for a multivariate mixed-Poisson model via copulas
# =========================================================


#' Fitting a parametric multivariate mixed-Poisson model 
#' 
#' Fitting a parametric multivariate mixed-Poisson model via a Gaussian copula with a log-normal marginal distribution of random effects 
#' 
#' @param formula a formula object with an obect of the type \code{Surv} on the left side and the terms on the right side. Note that \code{+strata()} is not supported.
#' @param data a 'data.frame' which has variables in 'formula'
#' @param frailty a charhacter string specifying a group 
#' @param dist a vector of character strings specifying a distribution for the baseline for each type of events. If it is a character string, it applies to all types of events. Default is "weibull". Other options are "ev", "loglogistic", and "lognormal". See \code{phreg}.  
#' @param ntype the number of types of events. 
#' @param model logical value: if \code{TRUE}, the model frame is returned.
#' @param model.matrix logical value: if \code{TRUE}, the x matirx is returned.
#' @param control an object of class specifying control options created by \code{\link{controlList}}.
#' @param init a list specifying inital values of a vector of variances of random effects and dependence parameters. Default initial value for a variance of random effect is 1 and 0 for the dependence parameter.
#' @param eventnames a vector of character string speficying event names
#' @param parallel logical value: if \code{TRUE}, it supports parallel execution.
#' @param ncore the number of cores if \code{parallel = TRUE}. 

#' @return an object of class \code{parmfrailty} representing the fit. 
#' \item{beta.coefficients}{a list of each types of events of the estimated regressoin coefficients.}
#' \item{logalpha.coefficients}{a list of each types of events of the estimated log of baseline parameters.}
#' \item{sig2}{a vector of the estimated variances of random effects.}
#' \item{rho}{a vector of the estimated dependence parameters from a Gaussian copula function.}
#' \item{ktau}{a vector of the estimated Kendall's tau from a Gaussian copula function.}
#' \item{logLik}{a value of log likelihood at the final values of the parameters.}
#' \item{iter}{the number of iterations used.}
#' \item{conv}{an integer code for the convergence. See \code{convergence} in \code{optim}.}
#' \item{var}{a variance-covariance matrix of the estimates.}
#' \item{betavar}{a list of each type of variances of the regression coefficient estimates.}
#' \item{logalphavar}{a list of each type of variances of the log of baseline parameter estimates.}
#' \item{thetavar}{a list of each type of variacnes of the regression coefficients and the log of baseline parameter estimates.}
#' \item{sig2var}{a vector of the variances of random effects variances estimates.}
#' \item{rhovar}{a vector of the dependence parameter estimates.}
#' \item{varnames}{a list of event names and variable names.}
#' \item{control}{a list of control arguments used.}
#' \item{n}{the number of sample size.}
#' \item{nevent}{a vector of the number of each types of events.}
#' \item{neventtype}{the number of types of events.}
#' \item{nknot}{the number nodes used in Gaussian-quadrature.}
#' \item{margins}{a vector of character strings specyting marginal distributions of random effects.}
#' The object will also contain the following: dist, copula, model, call, optionally x, and model.
#


#' @import eha
#' @export
parmfrailty <- function(formula, data, frailty = NULL, dist = "weibull", margins = "lnorm", ntype, model = FALSE, 
    model.matrix = FALSE, control = list(), init = list(), eventnames = NULL, parallel = FALSE, 
    ncore = NULL, ...) {
    
    
    
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
    
    
    if(length(dist) != J)  dist <- rep(dist, J)
    for (i in 1:J) {
        if (is.null(dist[i])) 
            stop("parametric distribution for event", i, "is missing")
    }
    
    Call <- match.call()
    
    if (missing(formula) | missing(data)) 
        stop("Missing arguments")
    
    if(length(margins) != J) margins <- rep(margins, J)
    cmargins <- margins
    margins <- ifelse(margins == "lnorm", 1, 2)
    
    mf <- Y <- type <- X <- coef.names <- id <- offset <- offset.names <- list()
    
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
        coef.names[[i]] <- colnames(X[[i]])
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
    
    #----- Initial parameters ----------------------------------------#
    
    newmodel <- newformula <- list()
    
    
    if (is.null(init$sig2)) {
        oldsig2 <- rep(1, J)
    } else oldsig2 <- init$sig2
    
    
    tmpcopula <- asCall0("normalCopula", NULL)
    
    if (is.null(init$rho)) {
        oldrho <- rep(0, K)
    } else {
        oldrho <- init$rho
    }
    
    # fit phreg ------------------------
    
    mphreg <- list()
    oldbeta <- oldalpha <- oldlam <- oldLam <- oldhazards <- oldHazards <- list()
    shape <- scale <- list()
    
    for (i in 1:J) {
        mphreg[[i]] <- eha::phreg(formula[[i]], data = data[[i]], dist = dist[i])
        oldbeta[[i]] <- mphreg[[i]]$coefficients[1:ncol(X[[i]])]
        shape[[i]] <- exp(as.numeric(mphreg[[i]]$coef[substr(names(mphreg[[i]]$coef), 
            5, 9) == "shape"]))
        scale[[i]] <- exp(as.numeric(mphreg[[i]]$coef[substr(names(mphreg[[i]]$coef), 
            5, 9) == "scale"]))
        
        
        oldalpha[[i]] <- c(shape[[i]], scale[[i]])
    }
    
    par = c(sapply(1:J, function(i) c(oldbeta[[i]], log(oldalpha[[i]]), log(oldsig2[i]))), 
        oldrho)
    
    
    ml <- try(optim(mPL, par = par, J = J, Y = Y, X = X, offset = offset, dist = dist, 
        margins = margins, ww = ww, vv = vv, nsample = nsample, N = N, id = id, nev_id = nev_id, 
        com = com, parallel = parallel, ncore = control$ncore, control = list(maxit = control$maxiter)))
    
    if (class(ml) == "try-error") {
        result <- list
        conv <- 1
        message("Fail to optimization")
    } else {
        
        # Estimates -------------------------------
        ncov = sapply(1:J, function(i) ncol(X[[i]]))
        
        est <- ml$par
        
        start <- 1
        est_beta <- est_alpha <- list()
        est_sig2 <- rep(0, J)
        
        for (i in 1:J) {
            est_beta[[i]] <- est[start:(start + ncov[i] - 1)]
            est_alpha[[i]] <- exp(est[(start + ncov[i]):(start + ncov[i] + 1)])
            est_sig2[i] <- exp(est[(start + ncov[i] + 2)])
            start <- start + ncov[i] + 3
        }
        
        
        est_rho <- est[-c(1:(sum(ncov) + 2 * J + J))]
        
        
        if (any((est_rho - 2 * control$delta) <= -1) | any((est_rho + 2 * control$delta) >= 
            1)) {
            message("Dependence parameter", " for event ", eventnames[((est_rho - 
                2 * control$delta) <= -1) | ((est_rho + 2 * control$delta) >= 1)], 
                " reaches the boundary", collapse = " ")
            conv2 <- 2
            result <- list(fail = c(est_beta, est_alpha, est_sig2, est_rho))
        } else {
            
            
            newhazards <- newHazards <- newlam <- newLam <- list()
            # Expectation
            for (i in 1:J) {
                shape = est_alpha[[i]][1]
                scale = est_alpha[[i]][2]
                if (dist[i] == "weibull") {
                  newhazards[[i]] <- eha::hweibull(Y[[i]][, 2], shape = shape, scale = scale)
                  newHazards[[i]] <- eha::Hweibull(Y[[i]][, 2], shape = shape, scale = scale) - 
                    eha::Hweibull(Y[[i]][, 1], shape = shape, scale = scale)
                } else if (dist[i] == "loglogistic") {
                  newhazards[[i]] <- eha::hllogis(Y[[i]][, 2], shape = shape, scale = scale[[i]])
                  newHazards[[i]] <- eha::Hllogis(Y[[i]][, 2], shape = shape, scale = scale) - 
                    eha::Hllogis(Y[[i]][, 1], shape = shape, scale = scale)
                } else if (dist[i] == "lognormal") {
                  sdlog <- 1/shape
                  meanlog <- log(scale)
                  
                  newhazards[[i]] <- eha::hlnorm(Y[[i]][, 2], meanlog = meanlog, 
                    sdlog = sdlog)
                  newHazards[[i]] <- eha::Hlnorm(Y[[i]][, 2], meanlog = meanlog, 
                    sdlog = sdlog) - eha::Hlnorm(Y[[i]][, 1], meanlog = meanlog, 
                    sdlog = sdlog)
                } else if (dist[i] == "ev") {
                  newhazards[[i]] <- eha::hEV(Y[[i]][, 2], shape = shape, scale = scale)
                  newHazards[[i]] <- eha::HEV(Y[[i]][, 2], shape = shape, scale = scale) - 
                    eha::HEV(Y[[i]][, 1], shape = shape, scale = scale)
                  
                } else {
                  stop(paste(dist, " is not availalble"))
                }
                newlam[[i]] <- (newhazards[[i]] * exp(X[[i]] %*% est_beta[[i]] + 
                  offset[[i]]))^(Y[[i]][, 3])
                newlam[[i]] <- aggregate(newlam[[i]], list(id[[i]]), prod, drop = T, 
                  simplify = T)[, -1]
                newLam[[i]] <- drop(rowsum(newHazards[[i]] * exp(X[[i]] %*% est_beta[[i]] + 
                  offset[[i]]), id[[i]], reorder = FALSE))
            }
            
            V <- parVar(beta = est_beta, alpha = est_alpha, sig2 = est_sig2, rho = est_rho, 
                lam = newlam, Lam = newLam, formula = formula, data = data, frailty = frailty, 
                J = J, margins = margins, num = num, control = control, ww = ww, 
                vv = vv, nsample = nsample, N = N, com = com, nev_id = nev_id, eventnames = eventnames, 
                parallel = parallel, dist = dist)
            
            
            p1_beta = unlist(lapply(1:J, function(i) length(est_beta[[i]])))
            p1_alpha = unlist(lapply(1:J, function(i) length(est_alpha[[i]])))
            p1_theta = unlist(lapply(1:J, function(i) length(est_beta[[i]]) + length(est_alpha[[i]])))
            
            p1 = sum(p1_theta)
            p2 = length(est_sig2)
            p3 = length(c(est_rho))
            
            if (anyNA(V)) {
                result <- list(fail = list(beta.coefficients = est_beta, alpha.coefficients = est_alpha, 
                  sig2 = est_sig2, rho = est_rho))
                message("no standard error is availble")
            } else {
                thetaV <- list()
                sig2V <- rep(0, J)
                p1_start = 1
                p1_end = -1
                diagV = diag(V)
                
                for (i in 1:J) {
                  p1_end = p1_end + p1_theta[i] + 1
                  thetaV[[i]] <- diagV[(p1_start):(p1_end)]
                  sig2V[i] <- diagV[(p1_end + 1)]
                  p1_start = p1_start + p1_theta[i] + 1
                }
                sig2V <- sig2V
                rhoV <- diagV[(1 + p1 + p2):(p1 + p2 + p3)]
                
                varNames <- list()
                varNames$eventnames <- eventnames
                varNames$theta.names <- lapply(1:J, function(i) paste(paste("event", 
                  eventnames[i], names(est_beta[[i]]), sep = "."), paste("log(shape)", 
                  paste("log(scale)"))))
                
                varNames$sig2.names <- paste("frailty", margins, "sd", sep = ".")
                varNames$rho.names <- "normalCopula"
                
                newktau <- rep(0, K)
                
                for (i in 1:K) {
                  tmpcopula <- asCall0("normalCopula", est_rho[i])
                  newktau[i] <- tau(eval(tmpcopula))
                }
                
                result <- list()
                result$beta.coefficients <- est_beta
                result$logalpha.coefficients <- lapply(1:J, function(i) log(est_alpha[[i]]))
                result$theta.coefficients <- lapply(1:J, function(i) c(est_beta[[i]], 
                  log(est_alpha[[i]])))
                for(i in 1:J) names(result$theta.coefficients[[i]]) <- c(coef.names[[i]], "log(shape)", "log(scale)")
                result$sig2 <- est_sig2
                result$rho <- est_rho
                result$ktau <- newktau
                result$logLik <- -ml$value
                result$iter <- ml$counts
                result$conv <- ml$convergence
                result$var <- V
                result$betavar <- lapply(1:J, function(i) thetaV[[i]][1:ncol(X[[i]])])
                result$logalphavar <- lapply(1:J, function(i) thetaV[[i]][-c(1:ncol(X[[i]]))])
                result$thetavar <- thetaV
                result$sig2var <- sig2V
                result$rhovar <- rhoV
		result$varnames <- varNames
            	result$control <- control
            	result$n <- nsample
            	result$nevent <- unlist(lapply(1:J, function(i) sum(nev_id[[i]])))
            	result$neventtype <- J
            	result$nknot <- control$nkont
            	result$margins <- cmargins
                result$offset <- offset
                
                if (model.matrix) 
                  result$x <- X
                if (model) 
                  result$model <- mf
                result$call <- Call
                result$dist <- dist
                result$copula <- "normalCopula"
                
                
            }
        }
    }
    
    attr(result, "class") <- "parmfrailty"
    
    return(result)
    
}
