\name{BMpoisson_covariates}
\alias{BMpoisson_covariates}
\title{Perform estimation on blockmodels for poisson probability distribution
aith covariates}
\description{
    With the provided network and blockmodel type, estimate number of groups,
parameters and node membership, and impact vector of covariates}
\usage{BMpoisson_covariates(membership_type, adj, verbosity = 4, plotlevel = 1)}
\arguments{
    \item{membership_type}{The type of node membership, e.g. "SBM", "LBM"}
    \item{adj}{The wieghted adjacency matrix}
    \item{verbosity}{The verbosity level, 0 means quiet, by default 4}
    \item{plotlevel}{More level is higher more plot of found values (as ICL
criterion), each time there are computed, are done. 0 means no plot.}
}
\examples{
    ## generation of one SBM network
    Z<-matrix(1,20,1)\%x\%diag(4)
    L<-70*matrix(runif(16),4,4)
    M_in_expectation_without_covariates<-Z\%*\%L\%*\%t(Z)
    Y1 <- matrix(runif(80*80),80,80)
    Y2 <- matrix(runif(80*80),80,80)
    M_in_expectation<-M_in_expectation_without_covariates*exp(4.2*Y1-1.2*Y2)
    M<-matrix(
        rpois(
            length(as.vector(M_in_expectation)),
            as.vector(M_in_expectation))
        ,80,80)

    ## estimation
    my_model <- BMpoisson_covariates("SBM",M,list(Y1,Y2),plotlevel=0)
    my_model$estimate()
    which.max(my_model$ICL)
}
