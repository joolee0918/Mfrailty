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

#include "Integrandf.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

# define XINFMAX DBL_MAX


using namespace Rcpp;

// Pairwise Likelihood

// Calculate E[ ]
//[[Rcpp::export()]]
NumericMatrix EU_cal(NumericVector cumsum1, NumericVector cumsum2,
  NumericVector dc, NumericVector ww1, NumericVector ww2,
  NumericVector vv1, NumericVector vv2, int nsample, int N, IntegerVector n1, IntegerVector n2, bool parallel, int ncore){
 
    int i,j;

    double sum1, sum2, sum3;
    

  //  Rcout<<dc<<"\n";
  //  Rcout<<cumsum1<<"\n";
  //  Rcout<<cumsum2<<"\n";
#ifdef _OPENMP
    omp_set_num_threads(ncore);
    int tid, nid, istart, iend;
#endif

  NumericMatrix Ef(nsample, 2);

#ifdef _OPENMP
#pragma omp parallel private (sum1, tid, i,j, istart, iend) shared(w, v, n, haz, cumsum, dg, N, res)
    {
        tid = omp_get_thread_num();
        nid = omp_get_num_threads();
        istart = tid*N/nid;
        iend = (tid+1)*N/nid;

  for ( i=0; i<nsample; i++){
      sum1=0.0, sum2=0.0, sum3=0.0;

      for(j=0; j<N; j++){
      sum1 +=  Integrand0(dc(j),
                          ww1(j), ww2(j), vv1(j), vv2(j), n1(i), n2(i), cumsum1(i), cumsum2(i));
      sum2 +=  Integrand1(dc(j),
                          ww1(j), ww2(j), vv1(j), vv2(j), n1(i), n2(i), cumsum1(i), cumsum2(i));
      sum3 +=  Integrand2(dc(j),
                          ww1(j), ww2(j), vv1(j), vv2(j), n1(i), n2(i), cumsum1(i), cumsum2(i));
    }
#pragma omp critical
      {
          Ef(i, 0) = sum2/sum1;
          Ef(i, 1) = sum3/sum1;
      }
    }
    }
    
#else
    for ( i=0; i<nsample; i++){
        sum1=0.0, sum2=0.0, sum3=0.0;

      for(j=0; j<N; j++){
        sum1 +=  Integrand0(dc(j),
                            ww1(j), ww2(j), vv1(j), vv2(j), n1(i), n2(i), cumsum1(i), cumsum2(i));
        sum2 +=  Integrand1(dc(j),
                            ww1(j), ww2(j), vv1(j), vv2(j), n1(i), n2(i), cumsum1(i), cumsum2(i));
        sum3 +=  Integrand2(dc(j),
                            ww1(j), ww2(j), vv1(j), vv2(j), n1(i), n2(i), cumsum1(i), cumsum2(i));
      }
        Ef(i, 0) = sum2/sum1;
        Ef(i, 1) = sum3/sum1;
        
    }
#endif
  return(Ef);
}

// Calculate E[ ]
//[[Rcpp::export()]]
double Ef_cal(NumericVector f, NumericVector cumsum1, NumericVector cumsum2,
               NumericVector dc, NumericVector ww1, NumericVector ww2,
                        NumericVector vv1, NumericVector vv2, int nsample, int N, IntegerVector n1, IntegerVector n2, bool parallel, int ncore){
    
    int i,j;

    double sum1, sum2;

    //  Rcout<<cumsum1<<"\n";
    //  Rcout<<cumsum2<<"\n";
#ifdef _OPENMP
    omp_set_num_threads(ncore);
    int tid, nid, istart, iend;
#endif
    
    double Ef = 0.0;
    
#ifdef _OPENMP
#pragma omp parallel private (sum1, tid, i,j, istart, iend) shared(w, v, n, haz, cumsum, dg, N, res, Ef)
    {
        tid = omp_get_thread_num();
        nid = omp_get_num_threads();
        istart = tid*N/nid;
        iend = (tid+1)*N/nid;
        
        for ( i=0; i<nsample; i++){
            sum1=0.0, sum2=0.0;
            
            for(j=0; j<N; j++){
                sum1 +=  Integrand0(dc(j),
                                    ww1(j), ww2(j), vv1(j), vv2(j), n1(i), n2(i), cumsum1(i), cumsum2(i));
                sum2 +=  Integrand3(f(j)), dc(j),
                                    ww1(j), ww2(j), vv1(j), vv2(j), n1(i), n2(i), cumsum1(i), cumsum2(i));
              
            }
#pragma omp critical
            {
                Ef += sum2/sum1;
            }
        }
    }
    
#else
    for ( i=0; i<nsample; i++){
        sum1=0.0, sum2=0.0;
        
        for(j=0; j<N; j++){
            sum1 +=  Integrand0(dc(j),
                                ww1(j), ww2(j), vv1(j), vv2(j), n1(i), n2(i), cumsum1(i), cumsum2(i));
            sum2 +=  Integrand3(f(j), dc(j),
                                ww1(j), ww2(j), vv1(j), vv2(j), n1(i), n2(i), cumsum1(i), cumsum2(i));
        }
        Ef += sum2/sum1;
        
    }
#endif
    return(Ef);
}

// Calculate E[ ]
//[[Rcpp::export()]]
NumericVector Ef1f2_cal(NumericVector f1, NumericVector f2, NumericVector cumsum1, NumericVector cumsum2,
			NumericVector dc, NumericVector ww1, NumericVector ww2,
                     NumericVector vv1, NumericVector vv2, int nsample, int N, IntegerVector n1, IntegerVector n2, bool parallel, int ncore){
    
    int i,j;
   
    double sum1, sum2, sum3;
    
   
    //  Rcout<<cumsum1<<"\n";
    //  Rcout<<cumsum2<<"\n";
#ifdef _OPENMP
    omp_set_num_threads(ncore);
    int tid, nid, istart, iend;
#endif
    
    NumericVector Ef(2);
    
#ifdef _OPENMP
#pragma omp parallel private (sum1, tid, i,j, istart, iend) shared(w, v, n, haz, cumsum, dg, N, res)
    {
        tid = omp_get_thread_num();
        nid = omp_get_num_threads();
        istart = tid*N/nid;
        iend = (tid+1)*N/nid;
        
        for ( i=0; i<nsample; i++){
            sum1=0.0, sum2=0.0, sum3=0.0;
            
            for(j=0; j<N; j++){
                sum1 +=  Integrand0(dc(j),
                                    ww1(j), ww2(j), vv1(j), vv2(j), n1(i), n2(i), cumsum1(i), cumsum2(i));
                sum2 +=  Integrand3(f1(j)), dc(j),
                                    ww1(j), ww2(j), vv1(j), vv2(j), n1(i), n2(i), cumsum1(i), cumsum2(i));
                sum3 +=  Integrand3(f2(j), dc(j),
                                    ww1(j), ww2(j), vv1(j), vv2(j), n1(i), n2(i), cumsum1(i), cumsum2(i));
            }
#pragma omp critical
            {
                Ef(0) += sum2/sum1;
                Ef(1) += sum3/sum1;
            }
        }
    }
    
#else
    for ( i=0; i<nsample; i++){
        sum1=0.0, sum2=0.0, sum3=0.0;
        
        for(j=0; j<N; j++){
            sum1 +=  Integrand0(dc(j),
                                ww1(j), ww2(j), vv1(j), vv2(j), n1(i), n2(i), cumsum1(i), cumsum2(i));
            sum2 +=  Integrand3(f1(j), dc(j),
                                ww1(j), ww2(j), vv1(j), vv2(j), n1(i), n2(i), cumsum1(i), cumsum2(i));
            sum3 +=  Integrand3(f2(j), dc(j),
                                ww1(j), ww2(j), vv1(j), vv2(j), n1(i), n2(i), cumsum1(i), cumsum2(i));
        }
        Ef(0) += sum2/sum1;
        Ef(1) += sum3/sum1;
        
    }
#endif
    return(Ef);
}

// Calculate Observed Likelihood

//[[Rcpp::export()]]
NumericVector ObslogL(NumericVector haz1, NumericVector haz2, NumericVector cumsum1, NumericVector cumsum2,
		      NumericVector dc,
  NumericVector ww1, NumericVector ww2, NumericVector vv1, NumericVector vv2, int nsample, int N, IntegerVector n1, IntegerVector n2, bool parallel, int ncore){

    
    int i,j;
   
    double sum1;
    
       
#ifdef _OPENMP
    omp_set_num_threads(ncore);
    int tid, nid, istart, iend;
#endif

  NumericVector res(nsample);

#ifdef _OPENMP
#pragma omp parallel private (sum1, sum2, tid, istart, iend, i, j) shared(w, v, n, cumsum, dg, N, Ef)
    {
        tid = omp_get_thread_num();
        nid = omp_get_num_threads();
        istart = tid*N/nid;
        iend = (tid+1)*N/nid;

  for ( i=0; i<nsample; i++){
    sum1=0.0;
    for(j=0; j<N; j++){
      sum1 +=  Integrand4(dc(j),
                          ww1(j), ww2(j), vv1(j), vv2(j), n1(i), n2(i), haz1(i), haz2(i), cumsum1(i), cumsum2(i));
    }
#pragma omp critical
      {
          res(i) = log(sum1);
      }
    }
    }

#else
    for ( i=0; i<nsample; i++){
        sum1=0.0;
    for(j=0; j<N; j++){
      sum1 +=  Integrand4(dc(j),
                          ww1(j), ww2(j), vv1(j), vv2(j), n1(i), n2(i), haz1(i), haz2(i), cumsum1(i), cumsum2(i));
    }
    res(i) = log(sum1);
  }
#endif
  return(res);
}


// Two-stage Pairwise Likelihood

// Calculate E[ ]
//[[Rcpp::export()]]
NumericVector EU_cal_ts( NumericVector cumsum,
                        double sig2, int margins, NumericVector w, NumericVector v,
                        int nsample, int N, NumericVector n, bool parallel, int ncore){
    
    int i,j;
    double sum1, sum2;

    NumericVector dg = dmargin_vec(v, margins, sig2, N);
    
#ifdef _OPENMP
    omp_set_num_threads(ncore);
    int tid, nid, istart, iend;
#endif
    NumericVector Ef(nsample);
    
#ifdef _OPENMP
#pragma omp parallel private (sum1, sum2, tid, istart, iend, i, j) shared(w, v, n, cumsum, dg, N, Ef)
    {
        tid = omp_get_thread_num();
        nid = omp_get_num_threads();
        istart = tid*N/nid;
        iend = (tid+1)*N/nid;
        
        for ( i=0; i<nsample; i++){
            sum1=0.0, sum2=0.0;
            for ( j=0; j<N; j++){
                sum1 +=  Integrand0_ts(dg(j),
                                    w(j),  v(j),  n(i),  cumsum(i));
                sum2 +=  Integrand1_ts(dg(j),
                                    w(j),  v(j),  n(i),  cumsum(i));
            }
#pragma omp critical
            {
                Ef(i) = sum2/sum1;
            }
        }
    }
    
#else
    for ( i=0; i<nsample; i++){
        sum1=0.0, sum2=0.0;
        for ( j=0; j<N; j++){
            sum1 +=  Integrand0_ts(dg(j),
                                w(j),  v(j),  n(i),  cumsum(i));
            sum2 +=  Integrand1_ts(dg(j),
                                w(j),  v(j),  n(i),  cumsum(i));
        }
        Ef(i) = sum2/sum1;
    }
#endif
    
   
    return(Ef);
}


// Calculate E[ ]
//[[Rcpp::export()]]
double Ef_cal_ts(NumericVector f, NumericVector cumsum,
			 double sig2, int margins, NumericVector w, NumericVector v,
			 int nsample, int N, NumericVector n, bool parallel, int ncore){
    
  int i,j;
  double sum1, sum2;

  NumericVector dg = dmargin_vec(v, margins, sig2, N);
    
#ifdef _OPENMP
  omp_set_num_threads(ncore);
  int tid, nid, istart, iend;
#endif
  double Ef=0.0;
    
#ifdef _OPENMP
#pragma omp parallel private (sum1, sum2, tid, istart, iend, i, j) shared(w, v, n, cumsum, dg, N, Ef)
  {
    tid = omp_get_thread_num();
    nid = omp_get_num_threads();
    istart = tid*N/nid;
    iend = (tid+1)*N/nid;
        
    for ( i=0; i<nsample; i++){
      sum1=0.0, sum2=0.0;
      for ( j=0; j<N; j++){
	sum1 +=  Integrand0_ts(dg(j),
			       w(j),  v(j),  n(i),  cumsum(i));
	sum2 +=  Integrand4_ts(f(j), dg(j),
			       w(j),  v(j),  n(i),  cumsum(i));
      }
#pragma omp critical
      {
	Ef += sum2/sum1;
      }
    }
  }
    
#else
  for ( i=0; i<nsample; i++){
    sum1=0.0, sum2=0.0;
    for ( j=0; j<N; j++){
      sum1 +=  Integrand0_ts(dg(j),
			     w(j),  v(j),  n(i),  cumsum(i));
      sum2 +=  Integrand4_ts(f(j), dg(j),
			     w(j),  v(j),  n(i),  cumsum(i));
    }
    Ef += sum2/sum1;
  }
#endif
    
   
  return(Ef);
}


// Calculate Observed Likelihood
//[[Rcpp::export()]]
NumericVector ObslogL_ts(NumericVector haz,  NumericVector cumsum,
                         double sig2, int margins, NumericVector w,   NumericVector v,  int nsample, int N, NumericVector n, bool parallel, int ncore){
    
    
    int i,j;
    double sum1;
    NumericVector dg = dmargin_vec(v, margins, sig2, N);
    
#ifdef _OPENMP
    omp_set_num_threads(ncore);
    int tid, nid, istart, iend;
#endif
    
    NumericVector res(nsample);
    
#ifdef _OPENMP
#pragma omp parallel private (sum1, tid, i,j, istart, iend) shared(w, v, n, haz, cumsum, dg, N, res)
    {
        tid = omp_get_thread_num();
        nid = omp_get_num_threads();
        istart = tid*N/nid;
        iend = (tid+1)*N/nid;
        
        for ( i=0; i<nsample; i++){
            sum1=0.0;
            for(j=0; j<N; j++){
                sum1 += Integrand3_ts(dg(j),
                                   w(j), v(j), n(i), haz(i), cumsum(i));
            }
#pragma omp critical
            {
                res(i) = log(sum1);
            }
        }
    }
#else
    for ( i=0; i<nsample; i++){
        sum1=0.0;
        for(j=0; j<N; j++){
            sum1 += Integrand3_ts(dg(j),
                               w(j), v(j), n(i), haz(i), cumsum(i));
        }
        res(i) = log(sum1);
    }
#endif
    return(res);
}



