
scalar_model <- setRefClass("scalar_model",
    contains = "model",
    fields = list(
        adj = "matrix"
    ),
    methods = list(
        number_of_nodes = function() { dim(adj) },
        network_to_cc = function() { list(adjacency = adj) },
        split_membership_model = function(Q)
        {
            membership <- memberships[[Q]]
            prediction <- predictions[[Q]]$adjacency
            error <- adj - prediction

            if(membership_name == "SBM")
            {
                result <- list()
                for(q in 1:Q)
                {
                    sub_classif <- coordinates_split(
                        cbind(error, t(error)),
                        membership$Z[,q]
                        )
                    Z <- cbind(membership$Z,membership$Z[,q])
                    Z[,q] <- Z[,q]*sub_classif
                    Z[,Q+1] <- Z[,Q+1]*(1-sub_classif)
                    result <- c(result, list(
                            getRefClass(membership_name)(from_cc=list(Z=Z))
                        ))
                }
                return(result)
            }
        },
        data_number = function()
        {
            if(membership_name=="SBM")
            {
                return(dim(adj)[1]*(dim(adj)[1]-1))
            }
            else
            {
                return(dim(adj)[1]*(dim(adj)[2]))
            }
        },
        precompute = function()
        {
            if(membership_name == "SBM")
            {
                if(length(precomputed)>0)
                {
                    return()
                }
                else
                {
                    if(length(predictions)!=0)
                    {
                        cat("comuptation of eigen decomposition used for initalizations")
                        error <- adj - predictions[[1]]$adjacency
                        W<- error %*% t(error)
                        W<-1/(1+exp(-W/sd(W)))
                        D<- diag(1/sqrt(rowSums(W)))
                        L<- D %*% W %*% D

                        precomputed$eigen <<- eigen(L, symmetric=TRUE)
                        
                        cat("\n")
                    }
                }
            }
        },
        provide_init = function(Q)
        {
            .self$precompute()
            if(membership_name == "SBM")
            {
                return(
                    SBM(
                        classif=lsbmkmeans(
                            precomputed$eigen$vectors[,1:Q],
                            Q
                        )
                    )
                )
            }
        }
    )
)


