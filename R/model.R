
setRefClass("model",
    fields = list(
        model_name = "character",       # e.g. "bernoulli"
        membership_name = "character",  # e.g. "SBM"

        autosave = "character", # autosave filename

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
        plotting = "character",

        # profiling
        profiling = "numeric",
        profiling_active = "logical",
        profiling_t = "numeric"
    ),
    methods = list(
        postinit = function()
        {
            if(membership_name != 'SBM'
               &&
               membership_name != 'SBM_sym'
               &&
               membership_name != 'LBM'
               )
            {
                stop(paste('Membership',membership_name,'unknown. Are you drunk?'))
            }
            if(length(verbosity)==0)
            {
                verbosity<<-4
            }

            if(length(profiling_active)==0)
            {
                profiling_active<<-FALSE
            }

            if(length(autosave)==0)
            {
                autosave<<-""
            }
        },
        save_now = function()
        {
            if(nzchar(autosave)>0)
            {
                saveRDS(.self,file=autosave)
            }
        },
        tic = function()
        {
            if(profiling_active)
            {
                profiling_t<<-cumtime()
            }
        },
        toc = function(field)
        {
            if(profiling_active)
            {
                t2 <- cumtime()
                if(is.na(profiling[field]))
                {
                    profiling[field] <<- 0
                }
                profiling[field] <<- profiling[field] + t2 - profiling_t
                profiling_t <<-t2
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
            if(length(plotting)==0)
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
                        list(list()),
                        2,reinitialization_effort)
                }
                else
                {
                    say(1,"Estimation for 1 groups")
                    do_with_inits(
                        list(getRefClass(membership_name)(
                            network_size=.self$number_of_nodes())),
                        list(list()),
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
                Qmin <- 2
            }
            else
            {
                Qmin <- 1
            }
            Q <- Qmin
            Qmax <- length(ICL)
            ret<-FALSE
            while(which.max(ICL)*1.5>length(ICL) || Q<Qmax || Q<Qmin+2)
            {
                Q<-Q+1

                say(3,"For",Q,"groups")

                if(Q>length(ICL) || changing_effort)
                {
                    say(4,"Init from spectral clustering")

                    tic()
                     
                    memberships_inits <- .self$provide_init(Q)
                    model_inits <- lapply(memberships_inits,function(x){list()})

                    toc('init_SC')
                    
                }
                else
                {
                    memberships_inits <- list()
                    model_inits <- list(list())
                }

                say(4,"Init from splitting groups from",Q-1,"groups")
               
                tic() 
                
                splitted<-.self$split_membership(Q-1)
                memberships_inits<-c(memberships_inits,splitted$memberships)
                model_inits<-c(model_inits,splitted$models)

                toc('init_split')

                if(length(memberships_inits)>0)
                {
                    r<-.self$do_with_inits(memberships_inits,model_inits,Q,reinitialization_effort)
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
            if(length(ICL)<=Qmin+1)
            {
                return(FALSE)
            }

            for(Q in seq(length(ICL)-1,Qmin+1))
            {
                say(3,"For",Q,"groups")
                say(4,"Init from merging groups from",Q+1,"groups")

                tic()

                merged <- merge_membership(Q+1)

                toc('init_merges')                

                if(length(merged$memberships)>0)
                {
                    r<-.self$do_with_inits(merged$memberships,merged$models,Q,reinitialization_effort)
                    ret<-ret||r
                }
                else
                {
                    say(4,"Already done")
                }

            }
            return(ret)
        },

        do_with_inits = function(memberships_inits,model_inits,Q,reinitialization_effort)
        {
            say(4,length(memberships_inits),"initializations provided")

            tic()
            
            filter<-sapply(
                memberships_inits,
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

            toc('estimation_already_tried')            

            nb_init_max <- floor(1+4*reinitialization_effort*sqrt(Q))

            say(5,length(memberships_inits)-sum(filter),"initializations already used")
            
            if(length(memberships_inits)>nb_init_max)
            {
                quality<-.self$membership_init_quality(memberships_inits,model_inits)
                seuil <- (-sort(-quality))[nb_init_max]
                filter <- filter & (quality >= seuil)
            }

            memberships_inits <- memberships_inits[filter]
            model_inits <- model_inits[filter]

            say(4,"Estimation with",length(memberships_inits),"initializations")

            ret<-FALSE

            
            if(length(memberships_inits)>0)
            {

                tic()

                results<-mcMap(
                    .self$do_one_estim,
                    memberships_inits,
                    model_inits,
                    mc.cores=detectCores(),
                    mc.preschedule=FALSE)
            
                toc('estimation_run')

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
               
                toc('estimation_computation_ICL')

                digest_already_tried <<- c(digest_already_tried,
                                    lapply(memberships_inits,function(x){x$digest()}))
                
                toc('estimation_adding_already_tried')

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

                    toc('estimation_saving_goods')

                    if(prod(dim(allICLs))==0)
                    {
                        allICLs<<-matrix(0,0,2)
                    }
                    allICLs<<-rbind(allICLs,c(Q,ICL[Q]))
                    
                    if(length(plotting)==0)
                    {
                        plot(allICLs)
                        points(ICL,type='b',col='red')
                    }
                    else
                    {
                        if(nzchar(plotting))
                        {
                            dev.pdf(plotting)
                            plot(allICLs)
                            points(ICL,type='b',col='red')
                            dev.off()
                        }
                    }



                
                }
                else
                {
                    say(5,"Useless, no better ICL criterion found")
                    say(5,"better ICL found:",max(ICLs))
                    say(5,"old ICL:",ICL[Q])
                }
            }

            .self$save_now()

            return(ret)
        },

        membership_init_quality = function(memberships_inits,model_inits)
        {
            tic()

            quals <- sapply(
                memberships_inits,
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
           
            toc('quality_already_computed') 
            
            
            if(any(is.na(quals)))
            {
                memberships_inits<-memberships_inits[is.na(quals)]
                model_inits<-model_inits[is.na(quals)]

                naquals<- simplify2array(
                    mcMap(
                        function(membership_init,model_init)
                        {
                            

                            r <- dispatcher(membership_name,
                                            membership_init$to_cc(),
                                            model_init,
                                            model_name,
                                            .self$network_to_cc(),
                                            FALSE)

                            local_ICL <- r$PL - .5*(r$model$n_parameters *
                                                log(.self$data_number())
                                        +
                                        getRefClass(membership_name)(
                                                from_cc=r$membership)$ICL_penalty())
                        },
                        memberships_inits,
                        model_inits
                    )
                )

                for(i in 1:length(memberships_inits))
                {
                    digest_already_quality_computed[[memberships_inits[[i]]$digest()]] <<- naquals[i]
                }

                quals[is.na(quals)] <- naquals
           
                toc('quality_computation') 
            }
            return(quals)
        },

        do_one_estim = function(membership_init,model_init)
        {
            return(
                dispatcher(membership_name,
                    membership_init$to_cc(),
                    model_init,
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
                splitted<-.self$split_membership_model(Q)
                digest_already_splitted <<- c(digest_already_splitted,list(d))
                return(splitted)
            }
        },
        merge_membership = function(Q)
        {
            membership <- memberships[[Q]]
            d<-membership$digest()
            if(any(sapply(digest_already_merged,function(x){x==d})))
            {
                return(list())
            }
            else
            {
                memberships_merged <- list()
                model_merged <- list()
                if(membership_name == "SBM" || membership_name == "SBM_sym")
                {
                    Z<-membership$Z
                    for(q in 1:(Q-1))
                    {
                        for(l in (q+1):Q)
                        {
                            Zn<-Z[,-l]
                            Zn[,q]<-Zn[,q]+Z[,l]
                            memberships_merged <- c(memberships_merged,list(
                                getRefClass(membership_name)(from_cc=list(Z=Zn))))
                            model_merged <- c(model_merged,list(
                                .self$merge_model(Q,q,l,sum(Z[,q]),sum(Z[,l]),0)))
                        }
                    }
                }
                else
                {
                    Q1<-ncol(membership$Z1)
                    Q2<-ncol(membership$Z2)
                    Z1 <- membership$Z1
                    Z2 <- membership$Z2
                    if(Q1>1)
                    {
                        for(q in 1:(Q1-1))
                        {
                            for(l in (q+1):Q1)
                            {
                                Zn<-as.matrix(Z1[,-l])
                                Zn[,q]<-Zn[,q]+Z1[,l]
                                memberships_merged <- c(memberships_merged,list(
                                    getRefClass(membership_name)(from_cc=list(Z1=Zn,Z2=Z2))))
                                model_merged <- c(model_merged,list(
                                    .self$merge_model(Q,q,l,sum(Z1[,q]),sum(Z1[,l]),1)))
                            }
                        }
                    }
                    if(Q2>1)
                    {
                        for(q in 1:(Q2-1))
                        {
                            for(l in (q+1):Q2)
                            {
                                Zn<-as.matrix(Z2[,-l])
                                Zn[,q]<-Zn[,q]+Z2[,l]
                                memberships_merged <- c(memberships_merged,list(
                                    getRefClass(membership_name)(from_cc=list(Z1=Z1,Z2=Zn))))
                                model_merged <- c(model_merged,list(
                                    .self$merge_model(Q,q,l,sum(Z2[,q]),sum(Z2[,l]),2)))
                            }
                        }
                    }
                }


                merged<-list(memberships=memberships_merged,models=model_merged)
                digest_already_merged <<- c(digest_already_merged,list(d))
                return(merged)
            }
        },
        precompute = function() {},
        plot_obs_pred = function(Q) {},
        plot_parameters = function(Q) {},
        split_model = function(Q,q,type) {list()},
        merge_model = function(Q,q,l,pq,pl,type) {list()},
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
            cat("blockmodels object\n")
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

