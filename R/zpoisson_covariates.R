
BM_poisson_covariates <- setRefClass("BM_poisson_covariates",
    contains = "scalar_model_with_covariates",
    methods = list(
        initialize = function(membership_type,adj,covariates,...)
        {
            if(class(covariates)=="matrix")
            {
                covariates_list <- list(covariates)
            }
            else
            {
                covariates_list <- covariates
            }
            .self$initFields(membership_name = membership_type,
                            model_name = "poisson_covariates",
                            adj = adj,
                            covariates = covariates_list,
                            ...)
            .self$postinit()
        },
        plot_transform = function(x) {log(1+x)},
        plot_parameters = function(Q)
        {
            matrixplot(.self$plot_transform(model_results[[Q]]$lambda))
        },
        prediction = function(Q)
        {
            B <- matrix(0,nrow(adj),ncol(adj))
            for(k in 1:length(covariates))
            {
                B <- B + model_results[[Q]]$beta[k] * covariates[[k]]
            }

            if(membership_name=='LBM')
            {
                return(
                    (
                        memberships[[Q]]$Z1
                        %*%
                        model_results[[Q]]$lambda
                        %*%
                        t(memberships[[Q]]$Z2)
                    ) * exp(B)
                )
            }
            else
            {
                return(
                    (
                        memberships[[Q]]$Z
                        %*%
                        model_results[[Q]]$lambda
                        %*%
                        t(memberships[[Q]]$Z)
                    ) * exp(B)
                )
            }
        },
        split_model=function(Q,q,type)
        {
            lambda <- model_results[[Q]]$lambda
            if(membership_name == "SBM" && membership_name == "SBM_sym")
            {
                lambda2<-matrix(1,Q+1,Q+1)
                lambda2[1:Q,1:Q] <- lambda
                lambda2[Q+1,]<-lambda[q,]
                lambda2[,Q+1]<-lambda[,q]
                lambda2[Q+1,Q+1]<-lambda[q,q]
            }
            else
            {
                Q1<-nrow(lambda)
                Q2<-ncol(lambda)
                if(type==1)
                {
                    lambda2<-matrix(1,Q1+1,Q2)
                    lambda2[1:Q1,]<-lambda
                    lambda2[Q1+1,]<-lambda[q,]
                }
                else
                {
                    lambda2<-matrix(1,Q1,Q2+1)
                    lambda2[,1:Q2]<-lambda
                    lambda2[,Q2+1]<-lambda[,q]
                }
            }
            return(list(lambda=lambda2,beta=model_results[[Q]]$beta))
        },
        merge_model=function(Q,q,l,pq,pl,type)
        {
            lambda<-model_results[[Q]]$lambda
            if(membership_name == "SBM" || membership_name == "SBM_sym")
            {
                lambda[q,] <- (pq*lambda[q,]+pl*lambda[l,])/(pq+pl)
                lambda[,q] <- (pq*lambda[,q]+pl*lambda[,l])/(pq+pl)
                lambda2<-lambda[-l,-l]
            }
            else
            {
                Q1<-nrow(lambda)
                Q2<-ncol(lambda)
                if(type==1)
                {
                    lambda[q,] <- (pq*lambda[q,]+pl*lambda[l,])/(pq+pl)
                    lambda2<-lambda[-l,]
                }
                else
                {
                    lambda[,q] <- (pq*lambda[,q]+pl*lambda[,l])/(pq+pl)
                    lambda2<-lambda[,-l]
                }
            }
            return(list(lambda=lambda2,beta=model_results[[Q]]$beta))
        }

    )

)

