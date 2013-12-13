
bernoulli <- setRefClass("bernoulli",
    contains = "scalar_model",
    methods = list(
        initialize = function(membership_type,adj,...)
        {
            .self$initFields(membership_name = membership_type,
                            model_name = "bernoulli",
                            adj = adj,
                            ...)
        },
        plot_parameters = function(Q)
        {
            matrixplot(model_results[[Q]]$pi)
        }
    )
)

