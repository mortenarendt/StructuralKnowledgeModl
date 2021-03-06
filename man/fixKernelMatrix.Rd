% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fixKernelMatrix.R
\name{fixKernelMatrix}
\alias{fixKernelMatrix}
\title{Modify similarity matrix to a smoothing kernel}
\usage{
fixKernelMatrix(P, thr = 1e-05, niter = 1000)
}
\arguments{
\item{P}{a n by n symmetric similarity matrix}

\item{thr}{convergency thresshold (default = 1e-5)}

\item{niter}{number of iterations (default = 1000)}
}
\description{
This functions calculates scales a similarity matrix P = P' such that P 1 = 1 and P'*1 = 1
}
