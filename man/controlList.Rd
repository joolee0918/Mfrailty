% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/controlList.R
\name{controlList}
\alias{controlList}
\title{Control parameters for emmfrailty and parmfrailty}
\usage{
controlList(control)
}
\arguments{
\item{tol.P}{Tolerance for convergence in the parameters. Default is 1e-04.}

\item{tol.L}{Tolerance for convergence in the log-likelihood. Default is 1e-10.}

\item{maxiter}{The maxiumum number of interations in Expectation-Maximization Algorithm for a semiparametric model and in \code{optim} for a parametric model.}

\item{delta}{maxiter An increament used for numerical differentiation in profile likelihood. Default is 1e-03.}

\item{nknot}{The number of nodes used in Gaussian-quadrature. Default is 20.}
}
\description{
Control parameters for emmfrailty and parmfrailty
}
