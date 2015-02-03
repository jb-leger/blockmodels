\name{BM_gaussian}
\alias{BM_gaussian}
\title{Perform estimation on blockmodels for gaussian probability distribution}
\description{
    With the provided network and blockmodel type, estimate number of groups,
parameters and node membership}
\usage{BM_gaussian(MEMBERSHIP_ARG, SCALAR_MODEL_ARG, MODEL_ARG)}
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
npc <- SBM_NPC # nodes per class
Q <- SBM_Q # classes
Z<-diag(Q)\%x\%matrix(1,npc,1)
Mu<-20*matrix(runif(Q*Q),Q,Q)
M<-matrix(rnorm(npc*Q*npc*Q,sd=10),npc*Q,npc*Q)+Z\%*\%Mu\%*\%t(Z) ## adjacency matrix

## estimation
my_model <- BM_gaussian("SBM",M)
my_model$estimate()
which.max(my_model$ICL)

##
## SBM symmetric
##

## generation of one SBM_sym network
npc <- SBM_NPC # nodes per class
Q <- SBM_Q # classes
Z<-diag(Q)\%x\%matrix(1,npc,1)
Mu<-20*matrix(runif(Q*Q),Q,Q)
Mu[lower.tri(Mu)]<-t(Mu)[lower.tri(Mu)]
M<-matrix(rnorm(npc*Q*npc*Q,sd=10),npc*Q,npc*Q)+Z\%*\%Mu\%*\%t(Z) ## adjacency matrix
M[lower.tri(M)]<-t(M)[lower.tri(M)]

## estimation
my_model <- BM_gaussian("SBM_sym",M)
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
Mu<-20*matrix(runif(Q[1]*Q[2]),Q[1],Q[2])
M<-matrix(rnorm(npc[1]*Q[1]*npc[2]*Q[2],sd=10),npc[1]*Q[1],npc[2]*Q[2])+Z1\%*\%Mu\%*\%t(Z2) ## adjacency matrix

## estimation
my_model <- BM_gaussian("LBM",M)
my_model$estimate()
which.max(my_model$ICL)
}
