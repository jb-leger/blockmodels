
model <- setRefClass("model",
    fields = list(
        model_name = "character",       # e.g. "bernoulli"
        membership_name = "character",  # e.g. "SBM"
        already_tried = "list",         # initialization already tried
        memberships = "list",           # found memberships
        model_results = "list",         # found model parameters
        PL = "numeric",                 # Pseudo liklihood of found models
        H = "numeric",                  # Entropy of found models
        ICL = "numeric",                # ICL of found models
        predictions = "list"            # Prediction of found models
    ),
    methods = list(
        do_estim = function()
        {
            if(length(memberships)==0)
            {
                do_with_inits(
                    list(getRefClass(membership_name)(
                        network_size=.self$number_of_nodes())),
                    1)
            }

            ICLmax = max(ICL)
            .self$estim_ascend()
            .self$estim_descend()

            while(TRUE)
            {
                if(max(ICL)==ICLmax)
                {
                    break
                }
                ICLmax = max(ICL)
                
                .self$estim_ascend()
                .self$estim_descend()
            }
        },

        estim_ascend = function()
        {
            Q <- 1
            Qmax <- length(ICL)
            while(which.max(ICL)*2>length(ICL) || Q<Qmax)
            {
                inits <- .self$split_membership(Q)
                Q<-Q+1
                if(length(inits)>0)
                {
                    .self$do_with_inits(inits,Q)
                }
                else
                {
                    break
                }
            }
        },

        estim_descend = function()
        {
            for(Q in seq(length(ICL)-1,2))
            {
                inits <- memberships[[Q+1]]$merges()

                if(length(inits)>0)
                {
                    .self$do_with_inits(inits,Q)
                }

            }
        },

        do_with_inits = function(inits,Q)
        {
            cat(paste("\rFor",Q,"groups with",
                      length(inits),"initalizations                 "))
            results<-mclapply(
                inits,
                .self$do_one_estim,
                mc.cores=detectCores())

            kmax <- which.max(sapply(results, function(x){x$PL}))
            r<-results[[kmax]]
            memberships[[Q]] <<-
                getRefClass(membership_name)(from_cc=r$membership)
            model_results[[Q]] <<- r$model
            predictions[[Q]] <<- r$prediction
            PL[Q] <<- r$PL
            H[Q] <<- r$H
            ICL[Q] <<- r$PL - .5*(r$model$n_parameters *
                log(.self$data_number()) + memberships[[Q]]$ICL_penalty())
            cat("\n\r                                                   \r")

        },

        do_one_estim = function(membership_init)
        {
            if(membership_init$into(already_tried))
            {
                return(FALSE)
            }
            dispatcher(membership_name,
                membership_init$to_cc(),
                model_name,
                .self$network_to_cc())
        }
    )
)

scalar_model <- setRefClass("scalar_model",
    contains = "model",
    fields = list(
        adj = "matrix"
    ),
    methods = list(
        number_of_nodes = function() { dim(adj) },
        network_to_cc = function() { list(adjacency = adj) },
        split_membership = function(Q)
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
        }

    )
)


