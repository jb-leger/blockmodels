
BM_poisson <- setRefClass("BM_poisson",
    contains = "scalar_model",
    methods = list(
        initialize = function(membership_type,adj,...)
        {
            .self$initFields(membership_name = membership_type,
                            model_name = "poisson",
                            adj = adj,
                            ...)
        },
        plot_parameters = function(Q)
        {
            matrixplot(model_results[[Q]]$lambda)
        },
        prediction = function(Q)
        {
            if(membership_name=='LBM')
            {
                return(
                    memberships[[Q]]$Z1
                     %*%
                    model_results[[Q]]$lambda
                     %*%
                    t(memberships[[Q]]$Z2)
                )
            }
            else
            {
                return(
                    memberships[[Q]]$Z
                     %*%
                    model_results[[Q]]$lambda
                     %*%
                    t(memberships[[Q]]$Z)
                )
            }
        }
    )
)

