
BM_bernoulli <- setRefClass("BM_bernoulli",
    contains = "scalar_model",
    methods = list(
        initialize = function(membership_type,adj,...)
        {
            .self$initFields(membership_name = membership_type,
                            model_name = "bernoulli",
                            adj = adj,
                            ...)
            .self$postinit()
        },
        plot_parameters = function(Q)
        {
            matrixplot(model_results[[Q]]$pi)
        },
        prediction = function(Q)
        {
            if(membership_name=='LBM')
            {
                return(
                    memberships[[Q]]$Z1
                     %*%
                    model_results[[Q]]$pi
                     %*%
                    t(memberships[[Q]]$Z2)
                )
            }
            else
            {
                return(
                    memberships[[Q]]$Z
                     %*%
                    model_results[[Q]]$pi
                     %*%
                    t(memberships[[Q]]$Z)
                )
            }
        }
    )
)

