
# membership definition

SBM <- setRefClass("SBM",
    fields = list(
        Z="matrix"
    ),
    methods = list(
        initialize = function(network_size=FALSE,classif=FALSE,from_cc=FALSE)
        {
            if(!classif[1])
            {
                if(network_size[1])
                {
                    Z <<- matrix(1,nrow=network_size[1],ncol=1)
                }
                else
                {
                    Z <<- from_cc$Z
                }
            }
            else
            {
                fclassif <- factor(classif)
                classif <- as.numeric(fclassif)
                Q <- length(levels(fclassif))
                Z <<- matrix(0,nrow=length(classif),ncol=Q)
                for(i in 1:length(classif))
                {
                    Z[i,classif[i]] <<- 1
                }
            }
        },
        into = function(memberships_list)
        {
            any(
                    sapply(
                        memberships_list,
                        function(x){
                            if(all(dim(x$Z)==dim(Z)))
                            {
                                return(all(x$Z == Z))
                            }
                            else
                            {
                                return(FALSE)
                            }
                        }
                    )
            )
        },
        into_v = function(membership_value_pairs_list)
        {
            bvec <- sapply(
                        membership_value_pairs_list,
                        function(x){
                            if(all(dim(x$membership$Z)==dim(Z)))
                            {
                                return(all(x$membership$Z == Z))
                            }
                            else
                            {
                                return(FALSE)
                            }
                        }

                    )
            if(any(bvec))
            {
                return(membership_value_pairs_list[[which(bvec)]]$ICL)
            }
            else
            {
                return(FALSE)
            }
        },
        to_cc = function()
        {
            list(Z=Z)
        },
        map = function()
        {
            list(
                C=apply(Z,1,which.max)
            )
        },
        ICL_penalty = function()
        {
            (dim(Z)[2]-1)*log(dim(Z)[1])
        },
        merges = function()
        {
            result <- list()
            Q <- dim(Z)[2]
            for(k1 in 1:(Q-1))
            {
                for(k2 in (k1+1):Q)
                {
                    Z2<-Z[,-k2]
                    Z2[,k1]<-Z[,k1]+Z[,k2]
                    result <- c(result,list(SBM(from_cc=list(Z=Z2))))
                }
            }
            return(result)
        }

    )
)
                

