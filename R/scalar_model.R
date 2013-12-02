
model <- setRefClass("model",
    fields = list(
        model_name = "character",       # e.g. "bernoulli"
        membership_name = "character",  # e.g. "SBM"
        already_tried = "list",         # initialization already tried
        already_quality_computed = "list", # initialization with quality already
                                           # computed
        already_splitted = "list",
        already_merged = "list",
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
                    1,'pass 0')
            }

            l<-TRUE
            n<-1
            while(l)
            {
            
                ra<-.self$estim_ascend(paste('pass',n))
                rb<-.self$estim_descend(paste('pass',n))
                l<-ra||rb
                n<-n+1
            }
        },

        estim_ascend = function(message)
        {
            Q <- 1
            Qmax <- length(ICL)
            ret<-FALSE
            while(which.max(ICL)*2>length(ICL) || Q<Qmax)
            {
                inits <- .self$split_membership(Q)
                Q<-Q+1
                if(length(inits)>0)
                {
                    r<-.self$do_with_inits(inits,Q,paste(message,'asc'))
                    ret<-ret||r
                }
                else
                {
                    if(Q>length(ICL))
                    {
                        break
                    }
                }
            }
            return(ret)
        },

        estim_descend = function(message)
        {
            ret<-FALSE
            for(Q in seq(length(ICL)-1,2))
            {
                inits <- merge_membership(memberships[[Q+1]])

                if(length(inits)>0)
                {
                    r<-.self$do_with_inits(inits,Q,paste(message,'desc'))
                    ret<-ret||r
                }

            }
            return(ret)
        },

        do_with_inits = function(inits,Q,message)
        {
            cat("\r")
            cat(paste(message,"for",Q,"groups with",
                      length(inits),"potential initalizations, "))

            filter <- sapply(inits,function(x){!x$into(.self$already_tried)})

            nb_init_max <- sqrt(16*Q)

            if(length(inits)>nb_init_max)
            {
                quality <- sapply(inits,.self$membership_init_quality)
                seuil <- (-sort(-quality))[floor(nb_init_max)]
                filter <- filter & (quality >= seuil)
            }

            inits <- inits[filter]
            cat(paste(length(inits),"are tried"))

            ret<-FALSE

            if(length(inits)>0)
            {

                results<-mclapply(
                    inits,
                    .self$do_one_estim,
                    mc.cores=detectCores())

                good <- FALSE

                ICLs <- sapply(results, function(r){
                            r$PL - .5*(r$model$n_parameters *
                                        log(.self$data_number())
                                +
                                getRefClass(membership_name)(
                                        from_cc=r$membership)$ICL_penalty())
                      })

                if(length(ICL)>=Q)
                {

                    if(ICL[Q]<max(ICLs))
                    {
                        good <- TRUE
                    }
                }
                else
                {
                    good <- TRUE
                }

                already_tried <<- c(already_tried, inits)

                if(good)
                {
                    kmax<-which.max(ICLs)

                    r<-results[[kmax]]
                    memberships[[Q]] <<-
                        getRefClass(membership_name)(from_cc=r$membership)
                    model_results[[Q]] <<- r$model
                    predictions[[Q]] <<- r$prediction
                    PL[Q] <<- r$PL
                    H[Q] <<- r$H
                    ICL[Q] <<- r$PL - .5*(r$model$n_parameters *
                        log(.self$data_number()) + memberships[[Q]]$ICL_penalty())
                    ret<-TRUE
                
                }
            }
            cat("\n")
            cat("\r                                                                                 \r")

            return(ret)
        },

        membership_init_quality = function(membership_init)
        {
            v<-membership_init$into_v(already_quality_computed)
            
            if(v)
            {
                return(v)
            }
            else
            {
                r <- dispatcher(membership_name,
                                membership_init$to_cc(),
                                model_name,
                                .self$network_to_cc(),
                                FALSE)

                local_ICL <- r$PL - .5*(r$model$n_parameters *
                                    log(.self$data_number())
                            +
                            getRefClass(membership_name)(
                                    from_cc=r$membership)$ICL_penalty())

                already_quality_computed<<-c(already_quality_computed,list(
                                        list(membership=membership_init,ICL=local_ICL)
                                        ))
                return(local_ICL)
            }
        },

        do_one_estim = function(membership_init)
        {
            return(
                dispatcher(membership_name,
                    membership_init$to_cc(),
                    model_name,
                    .self$network_to_cc(),
                    TRUE)
                )
        },
        split_membership = function(Q)
        {
            if(memberships[[Q]]$into(already_splitted))
            {
                return(list())
            }
            else
            {
                already_splitted <<- c(already_splitted, list(memberships[[Q]]))
                return(.self$split_membership_model(Q))
            }
        },
        merge_membership = function(membership)
        {
            if(membership$into(already_merged))
            {
                return(list())
            }
            else
            {
                already_merged <<- c(already_merged, list(membership))
                return(membership$merges())
            }
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
        }

    )
)


