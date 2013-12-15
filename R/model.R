
setRefClass("model",
    fields = list(
        model_name = "character",       # e.g. "bernoulli"
        membership_name = "character",  # e.g. "SBM"

        # digests are used here
        digest_already_tried = "list",     # initialization already tried
        digest_already_quality_computed = "list", 
                                           # initialization with quality already
                                           # computed each value is a list with
                                           # two members, first the digest,
                                           # second, the value
        digest_already_splitted = "list",
        digest_already_merged = "list",

        memberships = "list",           # found memberships
        model_results = "list",         # found model parameters
        PL = "numeric",                 # Pseudo liklihood of found models
        H = "numeric",                  # Entropy of found models
        ICL = "numeric",                # ICL of found models
        precomputed = "list",
        last_reinitialization_effort = "numeric",
        allICLs= "matrix",
        verbosity = "numeric",
        plotlevel = "numeric"
    ),
    methods = list(
        init_levels = function()
        {
            if(length(verbosity)==0)
            {
                verbosity<<-4
            }

            if(length(plotlevel)==0)
            {
                plotlevel<<-1
            }
        },
        say = function(level,...)
        {
            if(level<=verbosity)
            {
                if(level>1)
                {
                    for(i in 2:level)
                    {
                        cat("    ")
                    }
                }
                cat("-> ")
                cat(paste(...))
                cat("\n")
            }
        },
        estimate = function(reinitialization_effort=1)
        {
            .self$init_levels()
            
            if(plotlevel>=1)
            {
                dev.new()
            }
            if(!any(last_reinitialization_effort==reinitialization_effort))
            {
                digest_already_splitted <<- list()
                digest_already_merged <<- list()
                last_reinitialization_effort<<-reinitialization_effort
                changing_effort <- TRUE
            }
            else
            {
                changing_effort <- FALSE
            }
            
            if(length(memberships)==0)
            {
                if(membership_name=="LBM")
                {
                    say(1,"Estimation for 2 groups (1+1)")
                    do_with_inits(
                        list(getRefClass(membership_name)(
                            network_size=.self$number_of_nodes())),
                        2,reinitialization_effort)
                }
                else
                {
                    say(1,"Estimation for 1 groups")
                    do_with_inits(
                        list(getRefClass(membership_name)(
                            network_size=.self$number_of_nodes())),
                        1,reinitialization_effort)

                }
            }

            .self$precompute()

            l<-TRUE
            n<-1
            while(l)
            {
                say(1,"Pass",n)

                say(2,"With ascending number of groups")
                ra<-.self$estim_ascend(reinitialization_effort,changing_effort)

                say(2,"With descending number of groups")
                rb<-.self$estim_descend(reinitialization_effort)

                l<-ra||rb
                n<-n+1
                changing_effort<-FALSE
            }
        },

        estim_ascend = function(reinitialization_effort,changing_effort)
        {
            if(membership_name=="LBM")
            {
                Q <- 2
            }
            else
            {
                Q <- 1
            }
            Qmax <- length(ICL)
            ret<-FALSE
            while(which.max(ICL)*1.5>length(ICL) || Q<Qmax)
            {
                Q<-Q+1

                say(3,"For",Q,"groups")

                if(Q>length(ICL) || changing_effort)
                {
                    say(4,"Init from spectral clustering")
                    inits <- .self$provide_init(Q)
                }
                else
                {
                    inits <- list()
                }
                say(4,"Init from splitting groups from",Q-1,"groups")
                inits <- c(inits,.self$split_membership(Q-1))
                if(length(inits)>0)
                {
                    r<-.self$do_with_inits(inits,Q,reinitialization_effort)
                    ret<-ret||r
                }
                else
                {
                    say(4,"already done")
                    if(Q>length(ICL))
                    {
                        break
                    }
                }
            }
            return(ret)
        },

        estim_descend = function(reinitialization_effort)
        {
            ret<-FALSE
            if(membership_name=="LBM")
            {
                Qmin <- 2
            }
            else
            {
                Qmin <- 1
            }
            for(Q in seq(length(ICL)-1,Qmin+1))
            {
                say(3,"For",Q,"groups")
                say(4,"Init from merging groups from",Q+1,"groups")
                inits <- merge_membership(memberships[[Q+1]])

                if(length(inits)>0)
                {
                    r<-.self$do_with_inits(inits,Q,reinitialization_effort)
                    ret<-ret||r
                }
                else
                {
                    say(4,"Already done")
                }

            }
            return(ret)
        },

        do_with_inits = function(inits,Q,reinitialization_effort)
        {
            say(4,length(inits),"initializations provided")

            filter<-sapply(
                inits,
                function(init)
                {
                    d<-init$digest()
                    !any(
                        sapply(
                            digest_already_tried,
                            function(x)
                            {
                                x==d
                            }
                        )
                    )
                }
            )

            nb_init_max <- floor(1+4*reinitialization_effort*sqrt(Q))

            say(5,length(inits)-sum(filter),"initializations already used")
            
            if(length(inits)>nb_init_max)
            {
                quality<-.self$membership_init_quality(inits)
                seuil <- (-sort(-quality))[nb_init_max]
                filter <- filter & (quality >= seuil)
            }

            inits <- inits[filter]

            say(4,"Estimation with",length(inits),"initializations")

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

                digest_already_tried <<- c(digest_already_tried,
                                    lapply(inits,function(x){x$digest()}))


                if(good)
                {
                    say(5,"Better ICL criterion found")
                    say(5,"new ICL:",max(ICLs))
                    say(5,"old ICL:",ICL[Q])

                    kmax<-which.max(ICLs)

                    r<-results[[kmax]]
                    memberships[[Q]] <<-
                        getRefClass(membership_name)(from_cc=r$membership)

                    if(membership_name=="LBM")
                    {
                        say(5,memberships[[Q]]$show_short())
                    }

                    model_results[[Q]] <<- r$model
                    PL[Q] <<- r$PL
                    H[Q] <<- r$H
                    ICL[Q] <<- r$PL - .5*(r$model$n_parameters *
                        log(.self$data_number()) + memberships[[Q]]$ICL_penalty())
                    ret<-TRUE
                    if(prod(dim(allICLs))==0)
                    {
                        allICLs<<-matrix(0,0,2)
                    }
                    allICLs<<-rbind(allICLs,c(Q,ICL[Q]))
                    
                    if(plotlevel>=1)
                    {
                        plot(allICLs)
                        points(ICL,type='b',col='red')
                    }


                
                }
                else
                {
                    say(5,"Useless, no better ICL criterion found")
                    say(5,"better ICL found:",max(ICLs))
                    say(5,"old ICL:",ICL[Q])
                }
            }

            return(ret)
        },

        membership_init_quality = function(inits)
        {
            quals <- sapply(
                inits,
                function(init)
                {
                    qual <- digest_already_quality_computed[[init$digest()]]
                    if(is.null(qual))
                    {
                        return(NA)
                    }
                    else
                    {
                        return(qual)
                    }
                }
            )
            
            if(any(is.na(quals)))
            {
                inits<-inits[is.na(quals)]

                naquals<- sapply(
                    inits,
                    function(membership_init)
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
                    }
                )

                for(i in 1:length(inits))
                {
                    digest_already_quality_computed[[inits[[i]]$digest()]] <<- naquals[i]
                }

                quals[is.na(quals)] <- naquals
            }
            return(quals)
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
            d<-memberships[[Q]]$digest()
            if(any(sapply(digest_already_splitted,function(x){x==d})))
            {
                return(list())
            }
            else
            {
                splitted_membership<-.self$split_membership_model(Q)
                digest_already_splitted <<- c(digest_already_splitted,list(d))
                return(splitted_membership)
            }
        },
        merge_membership = function(membership)
        {
            d<-membership$digest()
            if(any(sapply(digest_already_merged,function(x){x==d})))
            {
                return(list())
            }
            else
            {
                merged_membership<-membership$merges()
                digest_already_merged <<- c(digest_already_merged,list(d))
                return(merged_membership)
            }
        },
        precompute = function() {},
        plot_obs_pred = function(Q) {},
        plot_parameters = function(Q) {},
        plot_all = function(Q = which.max(.self$ICL))
        {
            dev.new()
            memberships[[Q]]$plot()
            dev.new()
            .self$plot_obs_pred(Q)
            dev.new()
            .self$plot_parameters(Q)
        },
        show = function()
        {
            cat("lsbm object\n")
            cat(paste("    model:",model_name,"\n"))
            cat(paste("    membership:",membership_name,"\n"))
            cat(paste("    network:",.self$show_network(),"\n"))
            if(length(ICL)>0)
            {
                cat(paste("    maximum of ICL:",memberships[[which.max(ICL)]]$show_short(),"\n"))
                cat("    Most usefull fields and methods:\n")
                cat("        The following fields are indexed by the number of groups:\n")
                cat("            $ICL : vector of ICL\n")
                cat("            $PL : vector of pseudo log liklihood\n")
                cat("            $memberships : list of memberships founds by estimation\n")
                cat("            $model_results : models parameters founds by estimation\n")
                cat("        Estimation methods:\n")
                cat("            $estimate(reinitalization_effort=1) : to run again estimation with a\n")
                cat("                                                  higher reinitalization effort\n")
                cat("        Plotting methods:\n")
                cat("            $plot_obs_pred(Q) : to plot the obeserved and predicted network for Q groups\n")
                cat("            $plot_parameters(Q) : to plot the model_parameters for Q groups\n")
                cat("            $plot_all(Q=which.max(ICL)) : to plot memberships, and two previous plots\n")
            }
            else
            {
                cat("    Estimation not done.\n")
                cat("    Run $estimate(). You can specify a reinitialization effort, by default 1.\n")
            }

        }
    )
)

