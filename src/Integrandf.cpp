#define UMAX  1-1e-16
#define UMIN  1e-16
#define XEPS 1e-6

#include <Rcpp.h>

#define _USE_MATH_DEFINES
#include <math.h>
#include <float.h>
#include <Rmath.h>
#include "Rinternals.h"
#include <iostream>
#include <vector>


#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

# define XINFMAX DBL_MAX


using namespace Rcpp;


// Copula Function From VineCopula R packages
// Ulf Schepsmeier, Jakob Stoeber, Eike Christian Brechmann, Benedikt Graeler, Thomas Nagler and Tobias
// Erhardt (2018). VineCopula: Statistical Inference of Vine Copulas. R package version 2.1.8.
// https://CRAN.R-project.org/package=VineCopula


double StableGammaDivision(double x1, double x2)
{
    int i;
    double a1, a2, b1, b2, sum=1.0;
    a1 = fmod(MAX(x1,x2),1.0);
    a2 = MAX(x1,x2)-a1;
    b1 = fmod(MIN(x1,x2),1.0);
    b2 = MIN(x1,x2)-b1;
    if(a1==0.0 && b1==0.0)
    {
        for(i=1 ; i<(int)b2 ; i++) sum *= ((a1+a2)-(double)i)/((b1+b2)-(double)i);
        for(i=b2 ; i<(int)a2 ; i++) sum *= ((a1+a2)-(double)i);
    }
    else if(a1>0.0 && b1==0.0)
    {
        for(i=1 ; i<(int)b2 ; i++) sum *= ((a1+a2)-(double)i)/((b1+b2)-(double)i);
        for(i=(int)b2 ; i<=(int)a2 ; i++) sum *= ((a1+a2)-(double)i);
        sum *= gamma(a1);
    }
    else if(a1==0.0 && b1>0.0)
    {
        for(i=1 ; i<=(int)b2 ; i++) sum *= ((a1+a2)-(double)i)/((b1+b2)-(double)i);
        for(i=((int)b2+1) ; i<(int)a2 ; i++) sum *= ((a1+a2)-(double)i);
        sum /= gamma(b1);
    }
    else if(a1>0.0 && b1>0.0)
    {
        for(i=1 ; i<=(int)b2 ; i++) sum *= ((a1+a2)-(double)i)/((b1+b2)-(double)i);
        for(i=((int)b2+1) ; i<=(int)a2 ; i++) sum *= ((a1+a2)-(double)i);
        sum *= gamma(a1)/gamma(b1);
    }
    if(x2 > x1) sum = 1.0/sum;
    return sum;
}


// Vectorized version of log of copula density functions. Available copula: Independent, Gaussian, Student, Clayton, Gumbel, Frank
NumericVector logdCopula_vec(int family,  NumericVector u, NumericVector v, double theta, double nu, int n)
{
    int j;
    NumericVector res(n);
    double rho,  t1=0.0, t2=0.0, f;
    
    if(family==0) //independent
        res.fill(0);
    else if(family==1) //Gaussian
    {
        rho=theta;
        for(j=0;j<n;j++)
        {
            t1 = R::qnorm(u[j],0.0,1.0,1,0); t2 = R::qnorm(v[j],0.0,1.0,1,0);
            f = 1.0/sqrt(1.0-pow(rho,2.0))*exp((pow(t1,2.0)+pow(t2,2.0))/2.0+(2.0*rho*t1*t2-pow(t1,2.0)-pow(t2,2.0))/(2.0*(1.0-pow(rho,2.0))));
            if(log(f)>XINFMAX) res[j] = log(XINFMAX);
            else if(f < DBL_MIN) res[j] = log(DBL_MIN);
            else res[j] = log(f);
        }
    }
    else if(family==2) //Student
    {
        rho=theta;
        for(j=0;j<n;j++)
        {
            t1 = R::qt(u[j],nu,1,0); t2 = R::qt(v[j],nu,1,0);
            f = StableGammaDivision((nu+2.0)/2.0,nu/2.0)/(nu*M_PI*sqrt(1.0-pow(rho,2.0))*R::dt(t1,nu,0)*R::dt(t2,nu,0))*pow(1.0+(pow(t1,2.0)+pow(t2,2.0)-2.0*rho*t1*t2)/(nu*(1.0-pow(rho,2.0))),-(nu+2.0)/2.0);
            if(log(f)>XINFMAX) res[j] = log(XINFMAX);
            else if(f < DBL_MIN) res[j] = log(DBL_MIN);
            else res[j] = log(f);
        }
    }
    else if(family==3) //Clayton
    {
        if(theta == 0) res.fill(0);
        else if(theta < 1e-10) res.fill(0);
        else
        {
            for(j=0;j<n;j++)
            {
                f=log1p(theta)-(1.0+theta)*log(u[j]*v[j])-(2.0+1.0/(theta))*log(pow(u[j],-theta)+pow(v[j],-theta)-1.0);
                if(f>XINFMAX) res[j] = log(XINFMAX);
                else if(f<log(DBL_MIN)) res[j] = log(DBL_MIN);
                else res[j] = f;
            }
        }
    }
    else if(family==4) //Gumbel
    {
        for(j=0;j<n;j++)
        {
            t1 = pow(-log(u[j]),theta)+pow(-log(v[j]),theta);
            f= -pow(t1,1.0/(theta))+(2.0/(theta)-2.0)*log(t1)+(theta-1.0)*log(log(u[j])*log(v[j]))-log(u[j]*v[j])+log1p((theta-1.0)*pow(t1,-1.0/(theta)));
            
            if(f>XINFMAX) res[j] = log(XINFMAX);
            else if(f<log(DBL_MIN)) res[j] = log(DBL_MIN);
            else res[j] = f;
        }
    }
    else if(family==5) // Frank
    {
          if (fabs(theta) < 1e-10) {
              res.fill(0);
            } else {
                for(j=0;j<n;j++)
                {
                f = (theta*(exp(theta)-1.0)*exp(theta*v[j]+theta*u[j]+theta))/pow(exp(theta*v[j]+theta*u[j])-exp(theta*v[j]+theta)-exp(theta*u[j]+theta)+exp(theta),2.0);
                if(log(f)>XINFMAX) res[j] = log(XINFMAX);
                else if(f < DBL_MIN) res[j] = log(DBL_MIN);
                else res[j] = log(f);
            }
        }
    }
    
    return(res);
}


// copula density functions.
double logdCopula(int family, double u, double v, double theta, double nu)
{
    double res;
    double rho,  t1=0.0, t2=0.0, f;
    
       if(u<UMIN) u=UMIN;
        else if(u>UMAX) u=UMAX;
        if(v<UMIN) v=UMIN;
        else if(v>UMAX) v=UMAX;
  
    
    
    if(family==0) //independent
        res= 0;
    else if(family==1) //Gaussian
    {
        rho=theta;
        t1 = R::qnorm(u,0.0,1.0,1,0); t2 = R::qnorm(v,0.0,1.0,1,0);
            f = 1.0/sqrt(1.0-pow(rho,2.0))*exp((pow(t1,2.0)+pow(t2,2.0))/2.0+(2.0*rho*t1*t2-pow(t1,2.0)-pow(t2,2.0))/(2.0*(1.0-pow(rho,2.0))));
            if(log(f)>XINFMAX) res = log(XINFMAX);
            else if(f < DBL_MIN) res = log(DBL_MIN);
            else res = log(f);
        
    }
    else if(family==2) //Student
    {
        rho=theta;
            t1 = R::qt(u,nu,1,0); t2 = R::qt(v,nu,1,0);
            f = StableGammaDivision((nu+2.0)/2.0,nu/2.0)/(nu*M_PI*sqrt(1.0-pow(rho,2.0))*R::dt(t1,nu,0)*R::dt(t2,nu,0))*pow(1.0+(pow(t1,2.0)+pow(t2,2.0)-2.0*rho*t1*t2)/(nu*(1.0-pow(rho,2.0))),-(nu+2.0)/2.0);
            if(log(f)>XINFMAX) res = log(XINFMAX);
            else if(f < DBL_MIN) res = log(DBL_MIN);
            else res = log(f);
        
    }
    else if(family==3) //Clayton
    {
        if(theta == 0) res= 0;
        else if(theta < 1e-10) res= 0;
        else
        {
               f=log1p(theta)-(1.0+theta)*log(u*v)-(2.0+1.0/(theta))*log(pow(u,-theta)+pow(v,-theta)-1.0);
                if(f>XINFMAX) res = log(XINFMAX);
                else if(f<log(DBL_MIN)) res = log(DBL_MIN);
                else res = f;
            
        }
    }
    else if(family==4) //Gumbel
    {
           t1 = pow(-log(u),theta)+pow(-log(v),theta);
            f= -pow(t1,1.0/(theta))+(2.0/(theta)-2.0)*log(t1)+(theta-1.0)*log(log(u)*log(v))-log(u*v)+log1p((theta-1.0)*pow(t1,-1.0/(theta)));
            
            if(f>XINFMAX) res = log(XINFMAX);
            else if(f<log(DBL_MIN)) res = log(DBL_MIN);
            else res = f;
        
    }
    else  // Frank
    {
        if (fabs(theta) < 1e-10) {
            res = 0;
        } else {
               f = (theta*(exp(theta)-1.0)*exp(theta*v+theta*u+theta))/pow(exp(theta*v+theta*u)-exp(theta*v+theta)-exp(theta*u+theta)+exp(theta),2.0);
                if(log(f)>XINFMAX) res = log(XINFMAX);
                else if(f < DBL_MIN) res = log(DBL_MIN);
                else res = log(f);
         
        }
    }
    
    return(res);
}


// copula density functions.
NumericVector dCopula_vec(int family, NumericVector u, NumericVector v, double theta, double nu, int n)
{
    NumericVector res(n);
    double rho,  t1=0.0, t2=0.0, f;
    int j;

    if(family==0) //independent
      res.fill(1);
    else if(family==1) //Gaussian
    {
        rho=theta;
        if (theta == 0) {
            res.fill(1);
        } else {
            
        for(j=0;j<n;j++)
        {
        t1 = R::qnorm(u[j],0.0,1.0,1,0); t2 = R::qnorm(v[j],0.0,1.0,1,0);
        f = 1.0/sqrt(1.0-pow(rho,2.0))*exp((pow(t1,2.0)+pow(t2,2.0))/2.0+(2.0*rho*t1*t2-pow(t1,2.0)-pow(t2,2.0))/(2.0*(1.0-pow(rho,2.0))));
        res[j] = f;
        }
    }
    } else if(family==2) //Student
    {
        rho=theta;
        for(j=0;j<n;j++)
        {
        t1 = R::qt(u[j],nu,1,0); t2 = R::qt(v[j],nu,1,0);
        f = StableGammaDivision((nu+2.0)/2.0,nu/2.0)/(nu*M_PI*sqrt(1.0-pow(rho,2.0))*R::dt(t1,nu,0)*R::dt(t2,nu,0))*pow(1.0+(pow(t1,2.0)+pow(t2,2.0)-2.0*rho*t1*t2)/(nu*(1.0-pow(rho,2.0))),-(nu+2.0)/2.0);
        res[j] = f;
        }
    }
    else if(family==3) //Clayton
    {
        if(theta == 0) res.fill(1);
        else
        {
            for(j=0;j<n;j++)
            {
            f = (1.0+theta)*pow(u[j]*v[j],-1.0-theta)*pow(pow(u[j],-theta)+pow(v[j],-theta)-1.0,-2.0-1.0/(theta));
            f = MAX(f,0);
            res[j] = f;
            }
        }
    }
    else if(family==4) //Gumbel
    {
        if (theta == 1) {
            res.fill(1);
        } else {
            
        for(j=0;j<n;j++)
        {
        t1 = pow(-log(u[j]),theta)+pow(-log(v[j]),theta);
        t2 = exp(-pow(t1,1.0/(theta)));
        f = t2/(u[j]*v[j])*pow(t1,-2.0+2.0/(theta))*pow(log(u[j])*log(v[j]),theta-1.0)*(1.0+(theta-1.0)*pow(t1,-1.0/(theta)));
        res[j] = f;
        }
    }
    }else  // Frank
    {
        if (fabs(theta) < 1e-16) {
            res.fill(1);
        } else {
            for(j=0;j<n;j++)
            {
            f = (theta*(exp(theta)-1.0)*exp(theta*v[j]+theta*u[j]+theta))/pow(exp(theta*v[j]+theta*u[j])-exp(theta*v[j]+theta)-exp(theta*u[j]+theta)+exp(theta),2.0);
            res[j] = f;
            }
        }
    }
    
    return(res);
}



// Vectorized version of copula density functions.
double dCopula(int family, double u, double v, double theta, double nu)
{
    double res;
    double rho,  t1=0.0, t2=0.0, f;
    
    if(family==0) //independent
        res = 1;
    else if(family==1) //Gaussian
    {
        rho=theta;
        t1 = R::qnorm(u,0.0,1.0,1,0); t2 = R::qnorm(v,0.0,1.0,1,0);
        f = 1.0/sqrt(1.0-pow(rho,2.0))*exp((pow(t1,2.0)+pow(t2,2.0))/2.0+(2.0*rho*t1*t2-pow(t1,2.0)-pow(t2,2.0))/(2.0*(1.0-pow(rho,2.0))));
        res = f;
        
    }
    else if(family==2) //Student
    {
        rho=theta;
        t1 = R::qt(u,nu,1,0); t2 = R::qt(v,nu,1,0);
        f = StableGammaDivision((nu+2.0)/2.0,nu/2.0)/(nu*M_PI*sqrt(1.0-pow(rho,2.0))*R::dt(t1,nu,0)*R::dt(t2,nu,0))*pow(1.0+(pow(t1,2.0)+pow(t2,2.0)-2.0*rho*t1*t2)/(nu*(1.0-pow(rho,2.0))),-(nu+2.0)/2.0);
        res = f;
        
    }
    else if(family==3) //Clayton
    {
        if(theta == 0) res= 1;
        else if(theta < 1e-16) res= 1;
        else
        {
            f = (1.0+theta)*pow(u*v,-1.0-theta)*pow(pow(u,-theta)+pow(v,-theta)-1.0,-2.0-1.0/(theta));
            f = MAX(f,0);
            res = f;
            
        }
    }
    else if(family==4) //Gumbel
    {
        t1 = pow(-log(u),theta)+pow(-log(v),theta);
        t2 = exp(-pow(t1,1.0/(theta)));
        f = t2/(u*v)*pow(t1,-2.0+2.0/(theta))*pow(log(u)*log(v),theta-1.0)*(1.0+(theta-1.0)*pow(t1,-1.0/(theta)));
        res = f;
        
    }
    else  // Frank
    {
        if (fabs(theta) < 1e-16) {
            res = 0;
        } else {
            f = (theta*(exp(theta)-1.0)*exp(theta*v+theta*u+theta))/pow(exp(theta*v+theta*u)-exp(theta*v+theta)-exp(theta*u+theta)+exp(theta),2.0);
            res = f;
            
        }
    }
    
    return(res);
}


// Vectorized versinon of marginal distribution for random effect: log-normal, Gamma

double logdmargin (double u, int margin, double sig2){
    
    double logdG1;

if(margin==1) {
    logdG1 = R::dlnorm(u, -sig2/2, sqrt(sig2), 1.0);
  } else{
      logdG1 = R::dgamma(u, 1/sig2, sig2, 1.0);
 }
    return(logdG1);
}

double pmargin (double u, int margin, double sig2){
    
    double G1;
    if(margin==1) {
        G1 = R::plnorm(u, -sig2/2, sqrt(sig2), 1.0, 0.0);
    } else{
        G1 = R::pgamma(u, 1/sig2, sig2, 1.0, 0.0);
    }
    return(G1);
}


double dmargin (double u, int margin, double sig2){
    
    double dG1;
    
    if(margin==1) {
        dG1 = R::dlnorm(u, -sig2/2, sqrt(sig2), 0.0);
    } else{
        dG1 = R::dgamma(u, 1/sig2, sig2, 0.0);
    }
    return(dG1);
}


// marginal distribution for random effect: log-normal, Gamma


NumericVector pmargin_vec (NumericVector u, int margin, double sig2, int n){

    NumericVector G1(n);
    if(margin==1) {
        G1 = plnorm(u, -sig2/2, sqrt(sig2), 1.0, 0.0);
    } else{
        G1 = pgamma(u, 1/sig2, sig2, 1.0, 0.0);
    }
    return(G1);
}

NumericVector dmargin_vec (NumericVector u, int margin, double sig2, int n){
    
    NumericVector dG1(n);
     if(margin==1) {
        dG1 = dlnorm(u, -sig2/2, sqrt(sig2), 0.0);
    } else{
        dG1 = dgamma(u, 1/sig2, sig2, 0.0);
    }
    return(dG1);
}

// Integrand Functions for pairwise likelihood approach

// Denominator
double Integrand0(double dc_jk,
                  double w_j, double w_k, double u_j, double u_k, int n_j, int n_k, double cumLam_j, double cumLam_k) {

  double denom = w_j*w_k*pow(u_j, n_j)*exp((-1)*u_j*cumLam_j)*pow(u_k,n_k)*exp((-1)*u_k*cumLam_k)*dc_jk*exp(u_j+u_k);
  return(denom);
 }

// Numerator
/// EU1
 double Integrand1(double dc_jk,
                   double w_j, double w_k, double u_j, double u_k, int n_j, int n_k, double cumLam_j, double cumLam_k) {

  double u1num = w_j*w_k*u_j*pow(u_j,n_j)*exp((-1)*u_j*cumLam_j)*pow(u_k,n_k)*exp((-1)*u_k*cumLam_k)*dc_jk*exp(u_j+u_k);
  return(u1num);
   }

/// EU2
 double Integrand2(double dc_jk,
                   double w_j, double w_k, double u_j, double u_k, int n_j, int n_k, double cumLam_j, double cumLam_k) {

  double u2num = w_j*w_k*u_k*pow(u_j,n_j)*exp((-1)*u_j*cumLam_j)*pow(u_k,n_k)*exp((-1)*u_k*cumLam_k)*dc_jk*exp(u_j+u_k);
  return(u2num);
  }

/// E[value]
double Integrand3(double value, double dc_jk,
                  double w_j, double w_k, double u_j, double u_k, int n_j, int n_k, double cumLam_j, double cumLam_k) {

    double valuenum = w_j*w_k*value*pow(u_j,n_j)*exp((-1)*u_j*cumLam_j)*pow(u_k,n_k)*exp((-1)*u_k*cumLam_k)*dc_jk*exp(u_j+u_k);
 return(valuenum);
}

/// Observed Liklihood
double Integrand4(double dc_jk,
                  double w_j, double w_k, double u_j, double u_k, int n_j, int n_k, double haz_j, double haz_k, double cumLam_j, double cumLam_k) {
    
    double valuenum = w_j*w_k*pow(u_j,n_j)*haz_j*exp((-1)*u_j*cumLam_j)*pow(u_k,n_k)*haz_k*exp((-1)*u_k*cumLam_k)*dc_jk*exp(u_j+u_k);
 return(valuenum);
}


// Integrand Functions for two-stage pairwise likelihood approach

// Denominator
double Integrand0_ts(double dg, double w, double u, int n, double cumLam) {
    
    double denom = w*pow(u, n)*exp((-1)*u*cumLam)*dg*exp(u);
    return(denom);
}

// Numerator
/// EU
double Integrand1_ts(double dg, double w, double u, int n, double cumLam) {
    
    double u1num = w*u*pow(u, n)*exp((-1)*u*cumLam)*dg*exp(u);
    return(u1num);
}

/// E[value]
double Integrand2_ts(double newlogdg, double dg, double w, double u, int n, double cumLam){

    double valuenum =  w*newlogdg*pow(u, n)*exp((-1)*u*cumLam)*dg*exp(u);
    return(valuenum);
}

/// Observed Liklihood
double Integrand3_ts(double dg, double w, double u, int n, double haz, double cumLam){
    
    double valuenum =  w*haz*pow(u, n)*exp((-1)*u*cumLam)*dg*exp(u);
    return(valuenum);
}

double Integrand4_ts(double f, double dg, double w, double u, int n, double cumLam){

  double valuenum =  w*f*pow(u, n)*exp((-1)*u*cumLam)*dg*exp(u);
  return(valuenum);
}
