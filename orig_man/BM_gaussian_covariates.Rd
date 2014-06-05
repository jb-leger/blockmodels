\name{BM_gaussian_covariates}
\alias{BM_gaussian_covariates}
\title{Perform estimation on blockmodels for gaussian probability distribution
with covariates}
\description{
    With the provided network and blockmodel type, estimate number of groups,
parameters and node membership and impact vector of covariates}
\usage{BM_gaussian(MEMBERSHIP_ARG, SCALAR_MODEL_ARG, COVARIATES_ARG, MODEL_ARG)}
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
    Mu<-20*matrix(runif(16),4,4)
    Y1 <- matrix(runif(80*80),80,80)
    Y2 <- matrix(runif(80*80),80,80)
    M<-matrix(rnorm(80*80,sd=5),80,80)+Z\%*\%Mu\%*\%t(Z)+4.2*Y1-1.6*Y2 ## adjacency matrix

    ## estimation
    my_model <- BM_gaussian_covariates("SBM",M,list(Y1,Y2))
    my_model$estimate()
    which.max(my_model$ICL)

    ##
    ## SBM symmetric
    ##
    
    ## generation of one SBM_sym network
    Z<-matrix(1,20,1)\%x\%diag(4)
    Mu<-20*matrix(runif(16),4,4)
    Mu[lower.tri(Mu)]<-t(Mu)[lower.tri(Mu)]
    Y1 <- matrix(runif(80*80),80,80)
    Y2 <- matrix(runif(80*80),80,80)
    Y1[lower.tri(Y1)]<-t(Y1)[lower.tri(Y1)]
    Y2[lower.tri(Y2)]<-t(Y2)[lower.tri(Y2)]
    M<-matrix(rnorm(80*80,sd=5),80,80)+Z\%*\%Mu\%*\%t(Z)+4.2*Y1-1.6*Y2
    M[lower.tri(M)]<-t(M)[lower.tri(M)]

    ## estimation
    my_model <- BM_gaussian_covariates("SBM_sym",M,list(Y1,Y2))
    my_model$estimate()
    which.max(my_model$ICL)

    ##
    ## LBM
    ##

    ## generation of one LBM network
    Z1<-matrix(1,20,1)\%x\%diag(3)
    Z2<-matrix(1,25,1)\%x\%diag(6)
    Mu<-20*matrix(runif(18),3,6)
    Y1 <- matrix(runif(60*150),60,150)
    Y2 <- matrix(runif(60*150),60,150)
    M<-matrix(rnorm(60*150,sd=5),60,150)+Z1\%*\%Mu\%*\%t(Z2)+4.2*Y1-1.6*Y2 ## adjacency matrix

    ## estimation
    my_model <- BM_gaussian_covariates("LBM",M,list(Y1,Y2))
    my_model$estimate()
    which.max(my_model$ICL)
    
}
