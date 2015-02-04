require('blockmodels')
##
## SBM
##

## generation of one SBM network
npc <- 10 # nodes per class
Q <- 2 # classes
Z<-diag(Q)%x%matrix(1,npc,1)
Mu1<-8*matrix(runif(Q*Q),Q,Q)
Mu2<-8*matrix(runif(Q*Q),Q,Q)
M1<-matrix(rnorm(npc*Q*npc*Q,sd=5),npc*Q,npc*Q)+Z%*%Mu1%*%t(Z) ## adjacency
M2<-matrix(rnorm(npc*Q*npc*Q,sd=10),npc*Q,npc*Q)+Z%*%Mu2%*%t(Z) ## adjacency

## estimation
my_model <- BM_gaussian_multivariate_independent("SBM",list(M1,M2) , plotting='', explore_min=2, explore_max=2, ncores=2, verbosity=0)
my_model$estimate()
which.max(my_model$ICL)

##
## SBM symmetric
##

## generation of one SBM_sym network
npc <- 10 # nodes per class
Q <- 2 # classes
Z<-diag(Q)%x%matrix(1,npc,1)
Mu1<-8*matrix(runif(Q*Q),Q,Q)
Mu2<-8*matrix(runif(Q*Q),Q,Q)
Mu1[lower.tri(Mu1)]<-t(Mu1)[lower.tri(Mu1)]
Mu2[lower.tri(Mu2)]<-t(Mu2)[lower.tri(Mu2)]
M1<-matrix(rnorm(npc*Q*npc*Q,sd=5),npc*Q,npc*Q)+Z%*%Mu1%*%t(Z) ## adjacency
M2<-matrix(rnorm(npc*Q*npc*Q,sd=10),npc*Q,npc*Q)+Z%*%Mu2%*%t(Z) ## adjacency
M1[lower.tri(M1)]<-t(M1)[lower.tri(M1)]
M2[lower.tri(M2)]<-t(M2)[lower.tri(M2)]

## estimation
my_model <- BM_gaussian_multivariate_independent("SBM_sym",list(M1,M2) , plotting='', explore_min=2, explore_max=2, ncores=2, verbosity=0)
my_model$estimate()
which.max(my_model$ICL)

##
## LBM
##

## generation of one LBM network
npc <- c(20,10) # nodes per class
Q <- c(1,2) # classes
Z1<-diag(Q[1])%x%matrix(1,npc[1],1)
Z2<-diag(Q[2])%x%matrix(1,npc[2],1)
Mu1<-8*matrix(runif(Q[1]*Q[2]),Q[1],Q[2])
Mu2<-8*matrix(runif(Q[1]*Q[2]),Q[1],Q[2])
M1<-matrix(rnorm(npc[1]*Q[1]*npc[2]*Q[2],sd=5),npc[1]*Q[1],npc[2]*Q[2])+Z1%*%Mu1%*%t(Z2) ## adjacency
M2<-matrix(rnorm(npc[1]*Q[1]*npc[2]*Q[2],sd=10),npc[1]*Q[1],npc[2]*Q[2])+Z1%*%Mu2%*%t(Z2) ## adjacency

## estimation
my_model <- BM_gaussian_multivariate_independent("LBM",list(M1,M2) , plotting='', explore_min=2, explore_max=2, ncores=2, verbosity=0)
my_model$estimate()
which.max(my_model$ICL)
