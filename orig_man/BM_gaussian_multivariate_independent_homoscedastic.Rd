\name{BM_gaussian_multivariate_independant_homoscedastic}
\alias{BM_gaussian_multivariate_independant_homoscedastic}
\title{Perform estimation on blockmodels for multivariate independant
    homoscedastic gaussian probability distribution}
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
    Mu1<-matrix(c(2,0,0,2),2,2)\%x\%matrix(1,2,2)
    Mu2<-matrix(1,2,2)\%x\%matrix(c(2,0,0,2),2,2)
    M1<-matrix(rnorm(160*160,sd=5),160,160)+Z\%*\%Mu1\%*\%t(Z) ## adjacency
    M2<-matrix(rnorm(160*160,sd=5),160,160)+Z\%*\%Mu2\%*\%t(Z) ## adjacency

    ## estimation
    my_model <- BM_gaussian_multivariate_independent_homoscedastic("SBM",list(M1,M2))
    my_model$estimate()
    which.max(my_model$ICL)

    ##
    ## SBM symmetric
    ##
    
    ## generation of one SBM_sym network
    Z<-diag(4)\%x\%matrix(1,40,1)
    Mu1<-matrix(c(2,0,0,2),2,2)\%x\%matrix(1,2,2)
    Mu2<-matrix(1,2,2)\%x\%matrix(c(2,0,0,2),2,2)
    Mu1[lower.tri(Mu1)]<-t(Mu1)[lower.tri(Mu1)]
    Mu2[lower.tri(Mu2)]<-t(Mu2)[lower.tri(Mu2)]
    M1<-matrix(rnorm(160*160,sd=5),160,160)+Z\%*\%Mu1\%*\%t(Z) ## adjacency
    M2<-matrix(rnorm(160*160,sd=5),160,160)+Z\%*\%Mu2\%*\%t(Z) ## adjacency
    M1[lower.tri(M1)]<-t(M1)[lower.tri(M1)]
    M2[lower.tri(M2)]<-t(M2)[lower.tri(M2)]

    ## estimation
    my_model <- BM_gaussian_multivariate_independent_homoscedastic("SBM_sym",list(M1,M2))
    my_model$estimate()
    which.max(my_model$ICL)

    ##
    ## LBM
    ##

    ## generation of one LBM network
    Z1<-diag(2)\%x\%matrix(1,60,1)
    Z2<-diag(2)\%x\%matrix(1,120,1)
    Mu1<-matrix(c(2,0,2,0),2,2)
    Mu2<-matrix(c(2,2,0,0),2,2)
    M1<-matrix(rnorm(120*240,sd=5),120,240)+Z1\%*\%Mu1\%*\%t(Z2) ## adjacency
    M2<-matrix(rnorm(120*240,sd=5),120,240)+Z1\%*\%Mu2\%*\%t(Z2) ## adjacency

    ## estimation
    my_model <- BM_gaussian_multivariate_independent_homoscedastic("LBM",list(M1,M2))
    my_model$estimate()
    which.max(my_model$ICL)
    
}
