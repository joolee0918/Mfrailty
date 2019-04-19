library(Mfrailty)
library(copula)
tol <- .Machine$double.eps

formula = c(Surv(estart, estop, estatus)~trt+W,
            Surv(estart, estop, estatus)~trt+W,
            Surv(estart, estop, estatus)~trt+W)


nsample = 300; J =3 ; copula = "normalCopula"
tau = 1 ; probC = 0.1;
lam = c(1.5,1.75,2)
alpha = c(1,1.25,1.5)
beta11 = log(0.8); beta21 = log(0.8) ; beta31 = log(0.8);
beta12 = log(1.1); beta22 = log(1.1) ; beta32 = log(1.1);
sig2 = c(0.4^2, 0.4^2, 0.4^2)
rho = c(0.25, 0.25, 0.25)

margins = c("lnorm", "lnorm", "lnorm")
paircop = c("normalCopula")


sim = 1
data <- NULL
set.seed(sim)
data<- generatedata.f(nsample=nsample, J=J, copula=copula, margins=margins, tau=tau, probC=probC, sig21=sig2[1], lam1=lam[1], alp1=alpha[1], beta11=beta11, beta12=beta12,
       		      sig22=sig2[2], lam2=lam[2], alp2=alpha[2], beta21=beta21, beta22=beta22, sig23=sig2[3], lam3=lam[3], alp3=alpha[3], beta31=beta31, beta32=beta32,
         	      rho=rho, meanw=5, dispstr="un")

## pariwise likelihood ##
paircop=1
m1 <- emmfrailty(formula=formula,data= data, frailty="id", twostage=F, paircop=paircop, ntype=J,  parallel=T, ncore=8)

print(summary(m1))

#all.equal( m1$logLik,  -13962.09, tolerance = (10 ^ 14) * tol, scale = 1)
#all.equal( as.numeric(m1$beta.coefficients[[1]]), c(-0.29009018,  0.06017256), tolerance = (10 ^ 8) * tol)
#all.equal( as.numeric(m1$beta.coefficients[[2]]), c(-0.36034546,  0.02110336), tolerance = (10 ^ 8) * tol)
#all.equal( as.numeric(m1$beta.coefficients[[3]]), c(-0.26616888,  -0.01774453), tolerance = (10 ^ 8) * tol)
#all.equal( m1$sig2, c(0.16986106, 0.16983854, 0.13891057), tolerance = (10 ^ 9) * tol)
#all.equal( m1$rho, c(0.26680485, 0.38502940, 0.52134294), tolerance = (10 ^ 9) * tol)



## two-stage pairwise liklihood ##
paircop=1
m2 <- emmfrailty(formula=formula,data= data, frailty="id", twostage=T, paircop=paircop, ntype=J,  parallel=T, ncore=8)

print(summary(m2))

#all.equal( m2$logLik,  -13962.1, tolerance = (10 ^ 14) * tol, scale = 1)
#all.equal( as.numeric(m2$beta.coefficients[[1]]), c(-0.29027792, 0.06032894), tolerance = (10 ^ 8) * tol)
#all.equal( as.numeric(m2$beta.coefficients[[2]]), c(-0.36135267, 0.02116558), tolerance = (10 ^ 8) * tol)
#all.equal( as.numeric(m2$beta.coefficients[[3]]), c(-0.26515563, -0.01756548), tolerance = (10 ^ 8) * tol)
#all.equal( m2$sig2, c(0.16922662, 0.16921215, 0.13444554), tolerance = (10 ^ 9) * tol)
#all.equal( m2$rho, c(0.26774696, 0.38786461, 0.52380762), tolerance = (10 ^ 9) * tol)
