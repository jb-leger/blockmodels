\name{BM_gaussian_multivariate}
\alias{BM_gaussian_multivariate}
\title{Perform estimation on blockmodels for multivariate gaussian probability distribution}
\description{
    With the provided network and blockmodel type, estimate number of groups,
parameters and node membership}
\usage{BM_gaussian_multivariate(MEMBERSHIP_ARG, MULTIVARIATE_MODEL_ARG, MODEL_ARG)}
\arguments{
    MEMBERSHIP_ARG_ITEM
    MULTIVARIATE_MODEL_ARG_ITEM
    MODEL_ARG_ITEM
}
\examples{
npc <- SBM_NPC # nodes per class
Q <- SBM_Q # classes
Z<-diag(Q)\%x\%matrix(1,npc,1)
Mu1<-4*matrix(runif(Q*Q),Q,Q)
Mu2<-4*matrix(runif(Q*Q),Q,Q)
Noise1<-matrix(rnorm(npc*Q*npc*Q,sd=1),npc*Q,npc*Q)
Noise2<-matrix(rnorm(npc*Q*npc*Q,sd=1),npc*Q,npc*Q)
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
npc <- SBM_NPC # nodes per class
Q <- SBM_Q # classes
Z<-diag(Q)\%x\%matrix(1,npc,1)
Mu1<-4*matrix(runif(Q*Q),Q,Q)
Mu2<-4*matrix(runif(Q*Q),Q,Q)
Noise1<-matrix(rnorm(npc*Q*npc*Q,sd=1),npc*Q,npc*Q)
Noise2<-matrix(rnorm(npc*Q*npc*Q,sd=1),npc*Q,npc*Q)
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
npc <- LBM_NPC # nodes per class
Q <- LBM_Q # classes
Z1<-diag(Q[1])\%x\%matrix(1,npc[1],1)
Z2<-diag(Q[2])\%x\%matrix(1,npc[2],1)
Mu1<-4*matrix(runif(Q[1]*Q[2]),Q[1],Q[2])
Mu2<-4*matrix(runif(Q[1]*Q[2]),Q[1],Q[2])
Noise1<-matrix(rnorm(npc[1]*Q[1]*npc[2]*Q[2],sd=1),npc[1]*Q[1],npc[2]*Q[2])
Noise2<-matrix(rnorm(npc[1]*Q[1]*npc[2]*Q[2],sd=1),npc[1]*Q[1],npc[2]*Q[2])
M1<-Z1\%*\%Mu1\%*\%t(Z2) + Noise1 ## adjacency
M2<-Z1\%*\%Mu2\%*\%t(Z2) + 10*Noise1 + Noise2 ## adjacency

## estimation
my_model <- BM_gaussian_multivariate("LBM",list(M1,M2))
my_model$estimate()
which.max(my_model$ICL)
}
