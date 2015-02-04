##
## SBM
##

## generation of one SBM network
npc <- SBM_NPC # nodes per class
Q <- SBM_Q # classes
Z<-diag(Q)%x%matrix(1,npc,1)
L<-70*matrix(runif(Q*Q),Q,Q)
M_in_expectation_without_covariates<-Z%*%L%*%t(Z)
Y1 <- matrix(runif(npc*Q*npc*Q),npc*Q,npc*Q)
Y2 <- matrix(runif(npc*Q*npc*Q),npc*Q,npc*Q)
M_in_expectation<-M_in_expectation_without_covariates*exp(4.2*Y1-1.2*Y2)
M<-matrix(
    rpois(
        length(as.vector(M_in_expectation)),
        as.vector(M_in_expectation))
    ,npc*Q,npc*Q)

## estimation
my_model <- BM_poisson_covariates("SBM",M,list(Y1,Y2) EXAMPLE_TEST_ARGS)
my_model$estimate()
which.max(my_model$ICL)

##
## SBM symmetric
##

## generation of one SBM_sym network, we re-use one produced for SBM
npc <- SBM_NPC # nodes per class
Q <- SBM_Q # classes
Z<-diag(Q)%x%matrix(1,npc,1)
L<-70*matrix(runif(Q*Q),Q,Q)
L[lower.tri(L)]<-t(L)[lower.tri(L)]
M_in_expectation_without_covariates<-Z%*%L%*%t(Z)
Y1 <- matrix(runif(npc*Q*npc*Q),npc*Q,npc*Q)
Y2 <- matrix(runif(npc*Q*npc*Q),npc*Q,npc*Q)
Y1[lower.tri(Y1)]<-t(Y1)[lower.tri(Y1)]
Y2[lower.tri(Y2)]<-t(Y2)[lower.tri(Y2)]
M_in_expectation<-M_in_expectation_without_covariates*exp(4.2*Y1-1.2*Y2)
M<-matrix(
    rpois(
        length(as.vector(M_in_expectation)),
        as.vector(M_in_expectation))
    ,npc*Q,npc*Q)
M[lower.tri(M)]<-t(M)[lower.tri(M)]

## estimation
my_model <- BM_poisson_covariates("SBM_sym",M,list(Y1,Y2) EXAMPLE_TEST_ARGS)
my_model$estimate()
which.max(my_model$ICL)

##
## LBM
##

## generation of one LBM network
npc <- LBM_NPC # nodes per class
Q <- LBM_Q # classes
Z1<-diag(Q[1])%x%matrix(1,npc[1],1)
Z2<-diag(Q[2])%x%matrix(1,npc[2],1)
L<-70*matrix(runif(Q[1]*Q[2]),Q[1],Q[2])
M_in_expectation_without_covariates<-Z1%*%L%*%t(Z2)
Y1 <- matrix(runif(npc[1]*Q[1]*npc[2]*Q[2]),npc[1]*Q[1],npc[2]*Q[2])
Y2 <- matrix(runif(npc[1]*Q[1]*npc[2]*Q[2]),npc[1]*Q[1],npc[2]*Q[2])
M_in_expectation<-M_in_expectation_without_covariates*exp(4.2*Y1-1.2*Y2)
M<-matrix(
    rpois(
        length(as.vector(M_in_expectation)),
        as.vector(M_in_expectation))
    ,npc[1]*Q[1],npc[2]*Q[2])

## estimation
my_model <- BM_poisson_covariates("LBM",M,list(Y1,Y2) EXAMPLE_TEST_ARGS)
my_model$estimate()
which.max(my_model$ICL)
