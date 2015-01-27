\name{BM_gaussian_multivariate}
\alias{BM_gaussian_multivariate}
\title{Perform estimation on blockmodels for multivariate gaussian probability distribution}
\description{
    With the provided network and blockmodel type, estimate number of groups,
parameters and node membership}
\usage{BM_gaussian(MEMBERSHIP_ARG, MULTIVARIATE_MODEL_ARG, MODEL_ARG)}
\arguments{
    MEMBERSHIP_ARG_ITEM
    MULTIVARIATE_MODEL_ARG_ITEM
    MODEL_ARG_ITEM
}
\examples{

    ##
    ## SBM
    ##

    ## generation of one SBM network
    Z<-diag(4)\%x\%matrix(1,40,1)
    Mu1<-matrix(c(1,0,0,1),2,2)\%x\%matrix(1,2,2)
    Mu2<-matrix(1,2,2)\%x\%matrix(c(1,0,0,1),2,2)
    Noise1<-matrix(rnorm(160*160,sd=1),160,160)
    Noise2<-matrix(rnorm(160*160,sd=1),160,160)
    M1<- Z\%*\%Mu1\%*\%t(Z) + Noise1
    M2<- Z\%*\%Mu2\%*\%t(Z) + 10*Noise1 + Noise2

    ## estimation
    my_model <- BM_gaussian_multivariate("SBM",list(M1,M2))
    my_model$estimate()
    which.max(my_model$ICL)

    ##
    ## SBM symmetric
    ##
    
    ## generation of one SBM_sym network
    Z<-diag(4)\%x\%matrix(1,40,1)
    Mu1<-matrix(c(1,0,0,1),2,2)\%x\%matrix(1,2,2)
    Mu2<-matrix(1,2,2)\%x\%matrix(c(1,0,0,1),2,2)
    Noise1<-matrix(rnorm(160*160,sd=1),160,160)
    Noise2<-matrix(rnorm(160*160,sd=1),160,160)
    M1<- Z\%*\%Mu1\%*\%t(Z) + Noise1
    M2<- Z\%*\%Mu2\%*\%t(Z) + 10*Noise1 + Noise2
    M1[lower.tri(M1)]<-t(M1)[lower.tri(M1)]
    M2[lower.tri(M2)]<-t(M2)[lower.tri(M2)]

    ## estimation
    my_model <- BM_gaussian_multivariate("SBM_sym",list(M1,M2))
    my_model$estimate()
    which.max(my_model$ICL)

    ##
    ## LBM
    ##

    ## generation of one LBM network
    Z1<-diag(2)\%x\%matrix(1,60,1)
    Z2<-diag(2)\%x\%matrix(1,120,1)
    Mu1<-matrix(c(1,0,1,0),2,2)
    Mu2<-matrix(c(1,1,0,0),2,2)
    Noise1<-matrix(rnorm(120*240,sd=1),120,240)
    Noise2<-matrix(rnorm(120*240,sd=1),120,240)
    M1<-Z1\%*\%Mu1\%*\%t(Z2) + Noise1 ## adjacency
    M2<-Z1\%*\%Mu2\%*\%t(Z2) + 10*Noise1 + Noise2 ## adjacency

    ## estimation
    my_model <- BM_gaussian_multivariate("LBM",list(M1,M2))
    my_model$estimate()
    which.max(my_model$ICL)
    
}
