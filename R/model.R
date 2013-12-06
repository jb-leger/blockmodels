
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
        predictions = "list",           # Prediction of found models
        precomputed = "list",
        last_reinitialization_effort = "numeric"
    ),
    methods = list(
        do_estim = function(reinitialization_effort=1)
        {
            if(!any(last_reinitialization_effort==reinitialization_effort))
            {
                already_splitted <<- list()
                already_merged <<- list()
                last_reinitialization_effort<<-reinitialization_effort
                changing_effort <- TRUE
            }
            else
            {
                changing_effort <- FALSE
            }
            
            if(length(memberships)==0)
            {
                cat("for 1 group, ")
                do_with_inits(
                    list(getRefClass(membership_name)(
                        network_size=.self$number_of_nodes())),
                    1,reinitialization_effort)
                cat("\n")
            }

            .self$precompute()

            l<-TRUE
            n<-1
            while(l)
            {
            
                ra<-.self$estim_ascend(paste('pass',n),reinitialization_effort,changing_effort)
                rb<-.self$estim_descend(paste('pass',n),reinitialization_effort)
                l<-ra||rb
                n<-n+1
                changing_effort<-FALSE
            }
        },

        estim_ascend = function(message,reinitialization_effort,changing_effort)
        {
            Q <- 1
            Qmax <- length(ICL)
            ret<-FALSE
            while(which.max(ICL)*1.5>length(ICL) || Q<Qmax)
            {
                Q<-Q+1
                cat(paste(message,'asc, for',Q,'groups, '))
                if(Q>length(ICL) || changing_effort)
                {
                    inits <- list(.self$provide_init(Q))
                }
                else
                {
                    inits <- list()
                }
                inits <- c(inits,.self$split_membership(Q-1))
                if(length(inits)>0)
                {
                    r<-.self$do_with_inits(inits,Q,reinitialization_effort)
                    ret<-ret||r
                }
                else
                {
                    cat("already done")
                    if(Q>length(ICL))
                    {
                        break
                    }
                }
                cat("\n")
            }
            return(ret)
        },

        estim_descend = function(message,reinitialization_effort)
        {
            ret<-FALSE
            for(Q in seq(length(ICL)-1,2))
            {
                cat(paste(message,'desc, for',Q,'groups, '))
                inits <- merge_membership(memberships[[Q+1]])

                if(length(inits)>0)
                {
                    r<-.self$do_with_inits(inits,Q,reinitialization_effort)
                    ret<-ret||r
                }
                else
                {
                    cat("already done")
                }
                cat("\n")

            }
            return(ret)
        },

        do_with_inits = function(inits,Q,reinitialization_effort)
        {
            cat(paste("with",
                      length(inits),"initalizations, "))

            filter <- sapply(inits,function(x){!x$into(.self$already_tried)})

            nb_init_max <- floor(1+2*reinitialization_effort*sqrt(Q))

            cat(paste(sum(filter),"not already used, "))
            
            if(length(inits)>nb_init_max)
            {
                quality <- c(
                        mclapply(
                            inits,
                            .self$membership_init_quality,
                            mc.cores=detectCores(),
                            mc.preschedule=FALSE
                        ),
                        recursive=TRUE
                     )
                seuil <- (-sort(-quality))[nb_init_max]
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
                    mc.cores=detectCores(),
                    mc.preschedule=FALSE)

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
        },
        precompute = function() {}
    )
)

