##
## SBM
##

## generation of one SBM network
npc <- SBM_NPC # nodes per class
Q <- SBM_Q # classes
sigmo <- function(x){1/(1+exp(-x))}
Z<-diag(Q)%x%matrix(1,npc,1)
Mg<-8*matrix(runif(Q*Q),Q,Q)-4
Y1 <- matrix(runif(npc*Q*npc*Q),npc*Q,npc*Q)-.5
Y2 <- matrix(runif(npc*Q*npc*Q),npc*Q,npc*Q)-.5
M_in_expectation<-sigmo(Z%*%Mg%*%t(Z) + 5*Y1-3*Y2)
M<-1*(matrix(runif(npc*Q*npc*Q),npc*Q,npc*Q)<M_in_expectation)

## estimation
my_model <- BM_bernoulli_covariates("SBM",M,list(Y1,Y2) EXAMPLE_TEST_ARGS)
my_model$estimate()
which.max(my_model$ICL)


##
## SBM symmetric
##

## generation of one SBM_sym network
npc <- SBM_NPC # nodes per class
Q <- SBM_Q # classes
sigmo <- function(x){1/(1+exp(-x))}
Z<-diag(Q)%x%matrix(1,npc,1)
Mg<-8*matrix(runif(Q*Q),Q,Q)-4
Mg[lower.tri(Mg)]<-t(Mg)[lower.tri(Mg)]
Y1 <- matrix(runif(npc*Q*npc*Q),npc*Q,npc*Q)-.5
Y2 <- matrix(runif(npc*Q*npc*Q),npc*Q,npc*Q)-.5
Y1[lower.tri(Y1)]<-t(Y1)[lower.tri(Y1)]
Y2[lower.tri(Y2)]<-t(Y2)[lower.tri(Y2)]
M_in_expectation<-sigmo(Z%*%Mg%*%t(Z) + 5*Y1-3*Y2)
M<-1*(matrix(runif(npc*Q*npc*Q),npc*Q,npc*Q)<M_in_expectation)
M[lower.tri(M)]<-t(M)[lower.tri(M)]

## estimation
my_model <- BM_bernoulli_covariates("SBM_sym",M,list(Y1,Y2) EXAMPLE_TEST_ARGS)
my_model$estimate()
which.max(my_model$ICL)



##
## LBM
##

## generation of one LBM network
npc <- LBM_NPC # nodes per class
Q <- LBM_Q # classes
sigmo <- function(x){1/(1+exp(-x))}
Z1<-diag(Q[1])%x%matrix(1,npc[1],1)
Z2<-diag(Q[2])%x%matrix(1,npc[2],1)
Mg<-8*matrix(runif(Q[1]*Q[2]),Q[1],Q[2])-4
Y1 <- matrix(runif(npc[1]*Q[1]*npc[2]*Q[2]),npc[1]*Q[1],npc[2]*Q[2])-.5
Y2 <- matrix(runif(npc[1]*Q[1]*npc[2]*Q[2]),npc[1]*Q[1],npc[2]*Q[2])-.5
M_in_expectation<-sigmo(Z1%*%Mg%*%t(Z2) + 5*Y1-3*Y2)
M<-1*(matrix(runif(npc[1]*Q[1]*npc[2]*Q[2]),npc[1]*Q[1],npc[2]*Q[2])<M_in_expectation)

## estimation
my_model <- BM_bernoulli_covariates("LBM",M,list(Y1,Y2) EXAMPLE_TEST_ARGS)
my_model$estimate()
which.max(my_model$ICL)
