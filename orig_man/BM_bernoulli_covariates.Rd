\name{BM_bernoulli_covariates}
\alias{BM_bernoulli_covariates}
\title{Perform estimation on blockmodels for bernoulli probability distribution
aith covariates}
\description{
    With the provided network and blockmodel type, estimate number of groups,
parameters and node membership, and impact vector of covariates.}
\usage{BM_bernoulli_covariates(MEMBERSHIP_ARG, SCALAR_MODEL_ARG, COVARIATES_ARG, MODEL_ARG)}
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
         sigmo <- function(x){1/(1+exp(-x))}
         Z<-matrix(1,20,1)\%x\%diag(4)
         Mg<-4*matrix(runif(16),4,4)-2
         Y1 <- matrix(runif(80*80),80,80)-.5
         Y2 <- matrix(runif(80*80),80,80)-.5
         M_in_expectation<-sigmo(Z\%*\%Mg\%*\%t(Z) + .6*Y1-.4*Y2)
         M<-1*(matrix(runif(80*80),80,80)<M_in_expectation)
     
         ## estimation
         my_model <- BM_bernoulli_covariates("SBM",M,list(Y1,Y2))
         my_model$estimate()
         which.max(my_model$ICL)


         ##
         ## SBM symmetric
         ##
     
         ## generation of one SBM_sym network
         sigmo <- function(x){1/(1+exp(-x))}
         Z<-matrix(1,20,1)\%x\%diag(4)
         Mg<-4*matrix(runif(16),4,4)-2
         Mg[lower.tri(Mg)]<-t(Mg)[lower.tri(Mg)]
         Y1 <- matrix(runif(80*80),80,80)-.5
         Y2 <- matrix(runif(80*80),80,80)-.5
         Y1[lower.tri(Y1)]<-t(Y1)[lower.tri(Y1)]
         Y2[lower.tri(Y2)]<-t(Y2)[lower.tri(Y2)]
         M_in_expectation<-sigmo(Z\%*\%Mg\%*\%t(Z) + .6*Y1-.4*Y2)
         M<-1*(matrix(runif(80*80),80,80)<M_in_expectation)
         M[lower.tri(M)]<-t(M)[lower.tri(M)]
     
         ## estimation
         my_model <- BM_bernoulli_covariates("SBM_sym",M,list(Y1,Y2))
         my_model$estimate()
         which.max(my_model$ICL)



         ##
         ## LBM
         ##
     
         ## generation of one LBM network
         sigmo <- function(x){1/(1+exp(-x))}
         Z1<-matrix(1,40,1)\%x\%diag(2)
         Z2<-matrix(1,50,1)\%x\%diag(4)
         Mg<-4*matrix(runif(8),2,4)-2
         Y1 <- 2*matrix(runif(40*2*50*4),40*2,50*4)-1
         Y2 <- 2*matrix(runif(40*2*50*4),40*2,50*4)-1
         M_in_expectation<-sigmo(Z1\%*\%Mg\%*\%t(Z2) + .6*Y1-.4*Y2)
         M<-1*(matrix(runif(40*2*50*4),40*2,50*4)<M_in_expectation)
     
         ## estimation
         my_model <- BM_bernoulli_covariates("LBM",M,list(Y1,Y2))
         my_model$estimate()
         which.max(my_model$ICL)
}
