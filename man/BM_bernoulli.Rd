\name{BM_bernoulli}
\alias{BM_bernoulli}
\title{Perform estimation on blockmodels for bernoulli probability distribution}
\description{
    With the provided network and blockmodel type, estimate number of groups,
parameters and node membership}
\usage{BM_bernoulli(membership_type, adj, verbosity = 4, plotlevel = 1)}
\arguments{
    \item{membership_type}{The type of node membership, e.g. "SBM", "SBM_sym", "LBM"}
    \item{adj}{The adjacency matrix}
    \item{verbosity}{The verbosity level, 0 means quiet, by default 4}
    \item{plotlevel}{More level is higher more plot of found values (as ICL
criterion), each time there are computed, are done. 0 means no plot.}
}
\examples{
    ## generation of one SBM network
    Z<-matrix(1,20,1)\%x\%diag(4)
    P<-matrix(runif(16),4,4)
    M<-1*(matrix(runif(80*80),80,80)<Z\%*\%P\%*\%t(Z)) ## adjacency matrix

    ## estimation
    my_model <- BM_bernoulli("SBM",M,plotlevel=0)
    my_model$estimate()
    which.max(my_model$ICL)
}
