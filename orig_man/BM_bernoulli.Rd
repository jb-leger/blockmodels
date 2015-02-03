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
npc <- SBM_NPC # nodes per class
Q <- SBM_Q # classes
Z<-diag(Q)\%x\%matrix(1,npc,1)
P<-matrix(runif(Q*Q),Q,Q)
M<-1*(matrix(runif(npc*Q*npc*Q),npc*Q,npc*Q)<Z\%*\%P\%*\%t(Z)) ## adjacency matrix

## estimation
my_model <- BM_bernoulli("SBM",M)
my_model$estimate()
which.max(my_model$ICL)

##
## SBM symmetric
##

## generation of one SBM_sym network
npc <- SBM_NPC # nodes per class
Q <- SBM_Q # classes
Z<-diag(Q)\%x\%matrix(1,npc,1)
P<-matrix(runif(Q*Q),Q,Q)
P[lower.tri(P)]<-t(P)[lower.tri(P)]
M<-1*(matrix(runif(npc*Q*npc*Q),npc*Q,npc*Q)<Z\%*\%P\%*\%t(Z)) ## adjacency matrix
M[lower.tri(M)]<-t(M)[lower.tri(M)]

## estimation
my_model <- BM_bernoulli("SBM_sym",M)
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
P<-matrix(runif(Q[1]*Q[2]),Q[1],Q[2])
M<-1*(matrix(runif(npc[1]*Q[1]*npc[2]*Q[2]),npc[1]*Q[1],npc[2]*Q[2])<Z1\%*\%P\%*\%t(Z2)) ## adjacency matrix

## estimation
my_model <- BM_bernoulli("LBM",M)
my_model$estimate()
which.max(my_model$ICL)

}
