\name{BM_bernoulli}
\alias{BM_bernoulli}
\title{Perform estimation on blockmodels for bernoulli probability distribution}
\description{
    With the provided network and blockmodel type, estimate number of groups,
parameters and node membership}
\usage{BM_bernoulli(MEMBERSHIP_ARG, SCALAR_MODEL_ARG, MODEL_ARG)}
\arguments{
    MEMBERSHIP_ARG_ITEM
    SCALAR_MODEL_ARG_ITEM
    MODEL_ARG_ITEM
}
\examples{

    ##
    ## SBM
    ##

    ## generation of one SBM network
    Z<-matrix(1,20,1)\%x\%diag(4)
    P<-matrix(runif(16),4,4)
    M<-1*(matrix(runif(80*80),80,80)<Z\%*\%P\%*\%t(Z)) ## adjacency matrix

    ## estimation
    my_model <- BM_bernoulli("SBM",M)
    my_model$estimate()
    which.max(my_model$ICL)

    ##
    ## SBM symmetric
    ##
    
    ## generation of one SBM_sym network
    Z<-matrix(1,20,1)\%x\%diag(4)
    P<-matrix(runif(16),4,4)
    P[lower.tri(P)]<-t(P)[lower.tri(P)]
    M<-1*(matrix(runif(80*80),80,80)<Z\%*\%P\%*\%t(Z))
    M[lower.tri(M)]<-t(M)[lower.tri(M)]

    ## estimation
    my_model <- BM_bernoulli("SBM_sym",M)
    my_model$estimate()
    which.max(my_model$ICL)

    ##
    ## LBM
    ##

    ## generation of one LBM network
    Z1<-matrix(1,20,1)\%x\%diag(3)
    Z2<-matrix(1,25,1)\%x\%diag(6)
    P<-matrix(runif(18),3,6)
    M<-1*(matrix(runif(60*150),60,150)<Z1\%*\%P\%*\%t(Z2)) ## adjacency matrix

    ## estimation
    my_model <- BM_bernoulli("LBM",M)
    my_model$estimate()
    which.max(my_model$ICL)
    
}
