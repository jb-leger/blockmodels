\name{BM_poisson_covariates}
\alias{BM_poisson_covariates}
\title{Perform estimation on blockmodels for poisson probability distribution
aith covariates}
\description{
    With the provided network and blockmodel type, estimate number of groups,
parameters and node membership, and impact vector of covariates}
\usage{BM_poisson_covariates(MEMBERSHIP_ARG, SCALAR_MODEL_ARG, COVARIATES_ARG, MODEL_ARG)}
\arguments{
    MEMBERSHIP_ARG_ITEM
    SCALAR_MODEL_ARG_ITEM
    COVARIATES_ARG_ITEM
    MODEL_ARG_ITEM
}
\examples{

    ##
    ## SBM
    ##

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
    my_model <- BM_poisson_covariates("SBM",M,list(Y1,Y2))
    my_model$estimate()
    which.max(my_model$ICL)
    
    ##
    ## SBM symmetric
    ##

    ## generation of one SBM_sym network, we re-use one produced for SBM
    Z<-matrix(1,20,1)\%x\%diag(4)
    L<-70*matrix(runif(16),4,4)
    L[lower.tri(L)]<-t(L)[lower.tri(L)]
    M_in_expectation_without_covariates<-Z\%*\%L\%*\%t(Z)
    Y1 <- matrix(runif(80*80),80,80)
    Y2 <- matrix(runif(80*80),80,80)
    Y1[lower.tri(Y1)]<-t(Y1)[lower.tri(Y1)]
    Y2[lower.tri(Y2)]<-t(Y2)[lower.tri(Y2)]
    M_in_expectation<-M_in_expectation_without_covariates*exp(4.2*Y1-1.2*Y2)
    M<-matrix(
        rpois(
            length(as.vector(M_in_expectation)),
            as.vector(M_in_expectation))
        ,80,80)
    M[lower.tri(M)]<-t(M)[lower.tri(M)]

    ## estimation
    my_model <- BM_poisson_covariates("SBM_sym",M,list(Y1,Y2))
    my_model$estimate()
    which.max(my_model$ICL)

    ##
    ## LBM
    ##

    ## generation of one LBM network
    Z1<-matrix(1,20,1)\%x\%diag(3)
    Z2<-matrix(1,25,1)\%x\%diag(6)
    L<-70*matrix(runif(18),3,6)
    M_in_expectation_without_covariates<-Z1\%*\%L\%*\%t(Z2)
    Y1 <- matrix(runif(60*150),60,150)
    Y2 <- matrix(runif(60*150),60,150)
    M_in_expectation<-M_in_expectation_without_covariates*exp(4.2*Y1-1.2*Y2)
    M<-matrix(
        rpois(
            length(as.vector(M_in_expectation)),
            as.vector(M_in_expectation))
        ,60,150)

    ## estimation
    my_model <- BM_poisson_covariates("LBM",M,list(Y1,Y2))
    my_model$estimate()
    which.max(my_model$ICL)
 
}
