#ifndef INTEGRANDF_H
#define INTEGRANDF_H

#define UMAX  1-1e-16
#define UMIN  1e-16
#define XEPS 1e-6

#include <Rcpp.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

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

double StableGammaDivision(double x1, double x2);
NumericVector logdCopula_vec(int family,  NumericVector u, NumericVector v, double theta, double nu, int n);
double logdCopula(int family, double u, double v, double theta, double nu);
double logdmargin (double u, int margin, double sig2);
double dCopula(int family, double u, double v, double theta, double nu);
NumericVector dCopula_vec(int family, NumericVector u, NumericVector v, double theta, double nu, int n);
double pmargin (double u, int margin, double sig2);
double dmargin (double u, int margin, double sig2);
NumericVector pmargin_vec (NumericVector u, int margin, double sig2, int n);
NumericVector dmargin_vec (NumericVector u, int margin, double sig2, int n);
double Integrand0(double dc_jk,
                  double w_j, double w_k, double u_j, double u_k, int n_j, int n_k, double cumLam_j, double cumLam_k);
double Integrand1(double dc_jk,
                  double w_j, double w_k, double u_j, double u_k, int n_j, int n_k, double cumLam_j, double cumLam_k);
double Integrand2(double dc_jk,
                  double w_j, double w_k, double u_j, double u_k, int n_j, int n_k, double cumLam_j, double cumLam_k);
double Integrand3(double value, double dc_jk,
                  double w_j, double w_k, double u_j, double u_k, int n_j, int n_k, double cumLam_j, double cumLam_k);
double Integrand4(double dc_jk,
                  double w_j, double w_k, double u_j, double u_k, int n_j, int n_k, double haz_j, double haz_k, double cumLam_j, double cumLam_k);
double Integrand0_ts(double dg, double w, double u, int n, double cumLam);
double Integrand1_ts(double dg, double w, double u, int n, double cumLam);
double Integrand2_ts(double newlogdg, double dg, double w, double u, int n, double cumLam);
double Integrand3_ts(double dg, double w, double u, int n, double haz, double cumLam);
double Integrand4_ts(double f, double dg, double w, double u, int n, double cumLam);
#endif
