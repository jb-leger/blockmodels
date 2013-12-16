\name{BMpoisson}
\alias{BMpoisson}
\title{Perform estimation on blockmodels for poisson probability distribution}
\description{
    With the provided network and blockmodel type, estimate number of groups,
parameters and node membership}
\usage{BMbernoulli(membership_type, adj, verbosity = 4, plotlevel = 1)}
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
    M_in_expectation<-Z\%*\%L\%*\%t(Z)
    M<-matrix(
        rpois(
            length(as.vector(M_in_expectation)),
            as.vector(M_in_expectation))
        ,80,80)

    ## estimation
    my_model <- BMpoisson("SBM",M,plotlevel=0)
    my_model$estimate()
    which.max(my_model$ICL)
}
