\name{BM_gaussian}
\alias{BM_gaussian}
\title{Perform estimation on blockmodels for gaussian probability distribution}
\description{
    With the provided network and blockmodel type, estimate number of groups,
parameters and node membership}
\usage{BM_gaussian(membership_type, adj, verbosity=4, autosave="", plotting=character(0))}
\arguments{
    \item{membership_type}{The type of node membership, e.g. "SBM", "SBM_sym", "LBM"}
    \item{adj}{The adjacency matrix}
    \item{verbosity}{The verbosity level, 0 means quiet, by default 4}
    \item{autosave}{If \var{autosave} != "", after each estimation, the model
                    object is writed into file \var{autosave}. The model object
                    is readable by the function \var{readRDS}. Use-it for long
                    computation to allow restarting the estimation on system
                    crash.}
    \item{plotting}{Control plot of ICL values while the estimation is running.
        If plotting==character(0) (the default), plots are done on screen, if
            plotting=="", no plot are done, if plotting is a filename, plots are
    done in this filename.}
}
\examples{

    ##
    ## SBM
    ##

    ## generation of one SBM network
    Z<-matrix(1,40,1)\%x\%diag(4)
    Mu<-20*matrix(runif(16),4,4)
    M<-matrix(rnorm(160*160,sd=10),160,160)+Z\%*\%Mu\%*\%t(Z) ## adjacency matrix

    ## estimation
    my_model <- BM_gaussian("SBM",M)
    my_model$estimate()
    which.max(my_model$ICL)

    ##
    ## SBM symmetric
    ##
    
    ## generation of one SBM_sym network
    Z<-matrix(1,40,1)\%x\%diag(4)
    Mu<-20*matrix(runif(16),4,4)
    Mu[lower.tri(Mu)]<-t(Mu)[lower.tri(Mu)]
    M<-matrix(rnorm(160*160,sd=10),160,160)+Z\%*\%Mu\%*\%t(Z)
    M[lower.tri(M)]<-t(M)[lower.tri(M)]

    ## estimation
    my_model <- BM_gaussian("SBM_sym",M)
    my_model$estimate()
    which.max(my_model$ICL)

    ##
    ## LBM
    ##

    ## generation of one LBM network
    Z1<-matrix(1,40,1)\%x\%diag(3)
    Z2<-matrix(1,50,1)\%x\%diag(6)
    Mu<-20*matrix(runif(18),3,6)
    M<-matrix(rnorm(120*300,sd=20),120,300)+Z1\%*\%Mu\%*\%t(Z2) ## adjacency matrix

    ## estimation
    my_model <- BM_gaussian("LBM",M)
    my_model$estimate()
    which.max(my_model$ICL)
    
}
