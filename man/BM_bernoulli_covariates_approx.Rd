\name{BM_bernoulli_covariates_approx}
\alias{BM_bernoulli_covariates_approx}
\title{Perform estimation on blockmodels for bernoulli probability distribution
aith covariates}
\description{
    With the provided network and blockmodel type, estimate number of groups,
parameters and node membership, and impact vector of covariates. A parabolic
approximation is made to compute the liklihood, therefore the covariates effect
is supposed as not very high. In pratical usage this seems not be a problem.}
\usage{BM_bernoulli_covariates(membership_type, adj, covariates, verbosity=6, autosave='', plotting=character(0), exploration_factor=1.5, explore_min=4, explore_max=Inf, ncores=detectCores())}
\arguments{
    \item{membership_type}{The type of node membership, i.e. 'SBM', 'SBM_sym' or 'LBM'}
    \item{adj}{The adjacency matrix}
    \item{covariates}{Covariates matrix, or list of covariates matrices. Covariates matrix must have the same size than the adjacency matrix.}
    \item{verbosity}{The verbosity level, 0 means quiet. Level 1 display the phase of reinitialization. Level 2 display the level 1 and the ascending and descending phase for the number of groups. Level 3 display the level 2 and the number current number of groups which is estimated. Level 4 display the level 3 and the steps inside the estimation. Level 5 display the level 4, the current status of parallel running jobs and the current sub-step. Level 6 display level 5 and informations about ICL criteria found. Default is level 6. This parameter can be changed by accessing to the field $verbosity of the object.}
\item{autosave}{If \var{autosave} != '', after each estimation, the model object is writed into file \var{autosave}. The model object is readable by the function \var{readRDS}. Use-it for long computation to allow restarting the estimation on system crash. You can use it to alanyze the partial results when the estimation is running. This parameter can be changed by accessing to the field $autosave of the object.}
\item{plotting}{Control plot of ICL values while the estimation is running. If plotting==character(0) (the default), plots are done on screen, if plotting=='', no plot are done, if plotting is a filename, plots are done in this filename. This parameter can be changed by accessing the field $plotting of the object.}\item{exploration_factor}{Control the exploration of the number of groups. The exploration is stop when the number of groups reach exploration factor times the current maximum. By default 1.5. This parameter can be changed by accessing the field $exploration_factor of the object.}
\item{explore_min}{Explore to the explore_min number of groups even if the exploration_factor rule is satisfied. By default 4. This parameter can be changed by accessing the field $explore_min of the object.}
\item{explore_max}{Stop exploration after explore_max number of group in any case. By default Inf. This parameter can be changed by accessing the field $explore_max of the object.}
\item{ncores}{Number of parallel jobs to launch different EM intializations. By default detectCores(). This parameter can be changed by accessing the field $ncores of the object. This parameters is used only on Linux. Parallism is disabled on other plateform. (Not working on Windows, not tested on Mac OS, not tested on *BSD.)}

}
\examples{
         ##
         ## SBM
         ##
     
         ## generation of one SBM network         
         sigmo <- function(x){1/(1+exp(-x))}
         Z<-matrix(1,20,1)\%x\%diag(4)
         Mg<-8*matrix(runif(16),4,4)-4
         Y1 <- matrix(runif(80*80),80,80)-.5
         Y2 <- matrix(runif(80*80),80,80)-.5
         M_in_expectation<-sigmo(Z\%*\%Mg\%*\%t(Z) + .6*Y1-.4*Y2)
         M<-1*(matrix(runif(80*80),80,80)<M_in_expectation)
     
         ## estimation
         my_model <- BM_bernoulli_covariates_approx("SBM",M,list(Y1,Y2))
         my_model$estimate()
         which.max(my_model$ICL)


         ##
         ## SBM symmetric
         ##
     
         ## generation of one SBM_sym network
         sigmo <- function(x){1/(1+exp(-x))}
         Z<-matrix(1,20,1)\%x\%diag(4)
         Mg<-8*matrix(runif(16),4,4)-4
         Mg[lower.tri(Mg)]<-t(Mg)[lower.tri(Mg)]
         Y1 <- matrix(runif(80*80),80,80)-.5
         Y2 <- matrix(runif(80*80),80,80)-.5
         Y1[lower.tri(Y1)]<-t(Y1)[lower.tri(Y1)]
         Y2[lower.tri(Y2)]<-t(Y2)[lower.tri(Y2)]
         M_in_expectation<-sigmo(Z\%*\%Mg\%*\%t(Z) + .6*Y1-.4*Y2)
         M<-1*(matrix(runif(80*80),80,80)<M_in_expectation)
         M[lower.tri(M)]<-t(M)[lower.tri(M)]
     
         ## estimation
         my_model <- BM_bernoulli_covariates_approx("SBM_sym",M,list(Y1,Y2))
         my_model$estimate()
         which.max(my_model$ICL)



         ##
         ## LBM
         ##
     
         ## generation of one LBM network
         sigmo <- function(x){1/(1+exp(-x))}
         Z1<-matrix(1,40,1)\%x\%diag(2)
         Z2<-matrix(1,50,1)\%x\%diag(4)
         Mg<-8*matrix(runif(8),2,4)-4
         Y1 <- 2*matrix(runif(40*2*50*4),40*2,50*4)-1
         Y2 <- 2*matrix(runif(40*2*50*4),40*2,50*4)-1
         M_in_expectation<-sigmo(Z1\%*\%Mg\%*\%t(Z2) + .6*Y1-.4*Y2)
         M<-1*(matrix(runif(40*2*50*4),40*2,50*4)<M_in_expectation)
     
         ## estimation
         my_model <- BM_bernoulli_covariates_approx("LBM",M,list(Y1,Y2))
         my_model$estimate()
         which.max(my_model$ICL)
}
