##
## SBM
##

## generation of one SBM network
npc <- SBM_NPC # nodes per class
Q <- SBM_Q # classes
Z<-diag(Q)%x%matrix(1,npc,1)
Mu<-20*matrix(runif(Q*Q),Q,Q)
M<-matrix(rnorm(npc*Q*npc*Q,sd=10),npc*Q,npc*Q)+Z%*%Mu%*%t(Z) ## adjacency matrix

## estimation
my_model <- BM_gaussian("SBM",M EXAMPLE_TEST_ARGS)
my_model$estimate()
which.max(my_model$ICL)

##
## SBM symmetric
##

## generation of one SBM_sym network
npc <- SBM_NPC # nodes per class
Q <- SBM_Q # classes
Z<-diag(Q)%x%matrix(1,npc,1)
Mu<-20*matrix(runif(Q*Q),Q,Q)
Mu[lower.tri(Mu)]<-t(Mu)[lower.tri(Mu)]
M<-matrix(rnorm(npc*Q*npc*Q,sd=10),npc*Q,npc*Q)+Z%*%Mu%*%t(Z) ## adjacency matrix
M[lower.tri(M)]<-t(M)[lower.tri(M)]

## estimation
my_model <- BM_gaussian("SBM_sym",M EXAMPLE_TEST_ARGS)
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
Mu<-20*matrix(runif(Q[1]*Q[2]),Q[1],Q[2])
M<-matrix(rnorm(npc[1]*Q[1]*npc[2]*Q[2],sd=10),npc[1]*Q[1],npc[2]*Q[2])+Z1%*%Mu%*%t(Z2) ## adjacency matrix

## estimation
my_model <- BM_gaussian("LBM",M EXAMPLE_TEST_ARGS)
my_model$estimate()
which.max(my_model$ICL)
