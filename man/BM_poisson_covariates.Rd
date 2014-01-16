\name{BM_poisson_covariates}
\alias{BM_poisson_covariates}
\title{Perform estimation on blockmodels for poisson probability distribution
aith covariates}
\description{
    With the provided network and blockmodel type, estimate number of groups,
parameters and node membership, and impact vector of covariates}
\usage{BM_poisson_covariates(membership_type, adj, covariates, verbosity=4, autosave='', plotting=character(0))}
\arguments{
    \item{membership_type}{The type of node membership, i.e. 'SBM', 'SBM_sym' or 'LBM'}
    \item{adj}{The adjacency matrix}
    \item{covariates}{Covariates matrix, or list of covariates matrices. Covariates matrix must have the same size than the adjacency matrix.}
    \item{verbosity}{The verbosity level, 0 means quiet. Level 1 display the phase of reinitialization. Level 2 display the level 1 and the ascending and descending phase for the number of groups. Level 3 display the level 2 and the number current number of groups which is estimated. Level 4 display the level 3 and the steps inside the estimation. Level 5 display the level 4, the current status of parallel running jobs and the current sub-step. Level 6 display level 5 and informations about ICL criteria found. Default is level 6.}
\item{autosave}{If \var{autosave} != '', after each estimation, the model object is writed into file \var{autosave}. The model object is readable by the function \var{readRDS}. Use-it for long computation to allow restarting the estimation on system crash. You can use it to alanyze the partial results when the estimation is running.}
\item{plotting}{Control plot of ICL values while the estimation is running. If plotting==character(0) (the default), plots are done on screen, if plotting=='', no plot are done, if plotting is a filename, plots are done in this filename.}
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
