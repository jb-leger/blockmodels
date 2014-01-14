\name{BM_poisson}
\alias{BM_poisson}
\title{Perform estimation on blockmodels for poisson probability distribution}
\description{
    With the provided network and blockmodel type, estimate number of groups,
parameters and node membership}
\usage{BM_poisson(membership_type, SCALAR_MODEL_ARG, MODEL_ARG)}
\arguments{
    MEMBERSHIP_ARG_ITEM
    SCALAR_MODEL_ARG_ITEM
    MODEL_ARG_ITEM
}
\examples{
    
    #
    # SBM
    #

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
    my_model <- BM_poisson("SBM",M)
    my_model$estimate()
    which.max(my_model$ICL)

    ##
    ## SBM symmetric
    ##

    ## generation of one SBM_sym network
    Z<-matrix(1,20,1)\%x\%diag(4)
    L<-70*matrix(runif(16),4,4)
    L[lower.tri(L)]<-t(L)[lower.tri(L)]
    M_in_expectation<-Z\%*\%L\%*\%t(Z)
    M<-matrix(
        rpois(
            length(as.vector(M_in_expectation)),
            as.vector(M_in_expectation))
        ,80,80)
    M[lower.tri(M)]<-t(M)[lower.tri(M)]

    ## estimation
    my_model <- BM_poisson("SBM_sym",M)
    my_model$estimate()
    which.max(my_model$ICL)

    ##
    ## LBM
    ##

    ## generation of one LBM network
    Z1<-matrix(1,20,1)\%x\%diag(3)
    Z2<-matrix(1,25,1)\%x\%diag(6)
    L<-70*matrix(runif(18),3,6)
    M_in_expectation<-Z1\%*\%L\%*\%t(Z2)
    M<-matrix(
        rpois(
            length(as.vector(M_in_expectation)),
            as.vector(M_in_expectation))
        ,60,150)

    ## estimation
    my_model <- BM_poisson("LBM",M)
    my_model$estimate()
    which.max(my_model$ICL)
    
}
