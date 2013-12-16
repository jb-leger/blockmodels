
setRefClass("scalar_model_with_covariates",
    contains = "scalar_model",
    fields = list(
        covariates = "list"
    ),
    methods = list(
        network_to_cc = function()
        {
            list(
                adjacency = adj,
                covariates = covariates
            )
        }
    )
)
