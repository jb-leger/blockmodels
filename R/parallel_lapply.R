# parallel version of lapply, using non exported function of parallel package
# and compatible with the mclapply function of the parallel package, except the
# verbose arg.
#
# The difference with mclapply is:
# without prescedule
# with display
#
# In case of change of parallel package internal API, this function can be
# broken. A workaround is the following function (but we loose te display):
# parallel_lapply <- function(..., verbose=TRUE) 
# { 
#     return(parallel::mclapply(...))
# }
#

parallel_lapply <- function(X,FUN,...,mc.cores,mc.set.seed=FALSE,mc.silent=TRUE,verbose=TRUE)
{
    if (!is.vector(X) || is.object(X)) 
        X <- as.list(X)
    last_disp <- 0
    cat("\r")
    sx <- seq_along(X)
    res <- vector("list", length(sx))
        names(res) <- names(X)

    cores <- as.integer(mc.cores)

    if(.Platform$OS.type=="unix" && cores > 1)
    {
        do_parallel <- TRUE
    }
    else
    {
        do_parallel <- FALSE
    }

    if(do_parallel)
    {
        cores <- as.integer(mc.cores)

        if(length(X)<cores)
            cores <- length(X)


        ent <- rep(FALSE, length(X))
        fin <- rep(FALSE, length(X))
        jobid <- seq_len(cores)
        jobs <- lapply(jobid, function(i) mcparallel(FUN(X[[i]],...), mc.set.seed = mc.set.seed, silent = mc.silent))
        jobsp <- parallel:::processID(jobs)
        ent[jobid] <- TRUE
        has.errors <- 0L
        while (!all(fin)) {
            s <- parallel:::selectChildren(jobs, 0.5)
            if (is.null(s)) 
                break
            if (is.integer(s)) 
                for (ch in s) {
                    ji <- which(jobsp == ch)[1]
                    ci <- jobid[ji]
                    r <- parallel:::readChild(ch)
                    if (is.raw(r)) {
                        child.res <- unserialize(r)
                        if (inherits(child.res, "try-error")) 
                            has.errors <- has.errors + 1L
                        if (!is.null(child.res)) 
                            res[[ci]] <- child.res
                    }
                    else {
                        fin[ci] <- TRUE
                        if (!all(ent)) {
                            nexti <- which(!ent)[1]
                            jobid[ji] <- nexti
                            jobs[[ji]] <- mcparallel(FUN(X[[nexti]], 
                                                         ...), mc.set.seed = mc.set.seed, silent = mc.silent)
                            jobsp[ji] <- parallel:::processID(jobs[[ji]])
                            ent[nexti] <- TRUE
                        }
                    }
                }

            if(verbose)
            {
                cat("\r")
                for(theta in 1:last_disp)
                {
                    cat(" ")
                }
                nb_jobs_done <- sum(fin)
                nb_jobs_running <- min(cores, length(X)-nb_jobs_done)
                nb_jobs_pending <- length(X)-nb_jobs_running-nb_jobs_done
                cat("\r")
                message <- paste(
                    "Parallel jobs:",
                    nb_jobs_done,
                    "done,",
                    nb_jobs_running,
                    "running,",
                    nb_jobs_pending,
                    "pending "
                )
                last_disp <- nchar(message)
                cat(message)
            }

        }

        if (has.errors)
        {
            if(verbose)
                cat("\n")
            warning(gettextf("%d function calls resulted in an error", 
                             has.errors), domain = NA)
        }
        
    }
    else
    {
        k<-0
        nb_jobs_total <- length(sx)
        nb_jobs_done <- 0
        for(x in sx)
        {
            if(verbose)
            {
                cat("\r")
                for(theta in 1:last_disp)
                {
                    cat(" ")
                }
                cat("\r")
                message <- paste(
                    'Non-parallel jobs: ',
                    nb_jobs_done,
                    '/',
                    nb_jobs_total,
                    ' ',
                    sep=''
                )
                last_disp <- nchar(message)
                cat(message)
            }
            
            res[[x]] <- FUN(X[[x]],...)
            nb_jobs_done <- nb_jobs_done + 1
        }
    }
    
    if(verbose)
    {
        cat("\r")    
        for(theta in 1:last_disp)
        {
            cat(" ")
        }
        cat("\r")
    }


    return(res)
}



