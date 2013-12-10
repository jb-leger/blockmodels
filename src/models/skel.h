
/*
 * This is package skeleton for MODEL model
 *
 * You should delete this header when the MODEL model is ready for use to avoid
 * any confusion
 */

class MODEL
{
    public:

    class network
    {
        public:
        /* Here you should put all the variable which define the network
         * for scalar model, a mat field is good choice.
         * See bernoulli, poisson model for example.
         *
         * Covariates must also stocked here.
         * See poisson_covariates model for example
         */


        /* Here you should add all precomputed values which depends only on the
         * network usefull in various functions, to avoid computing these value
         * many time.
         */

        network(Rcpp::List & network_from_R)
        {
            /* Here you must define how initialize the fields describing the
             * network. Precomputed value must be computed here.
             * See bernoulli, poisson or poisson_covariates for example.
             *
             * For scalar network, the provided list have a adjacency field
             * which contains the adjacency matrix.
             * For scalar network with covariates vectors on edges (existent or
             * not), the provided list have a covariates field which contain a
             * list of matrices. This matrices have the same size than the
             * adjacency matrix, and the i-th matrix is the matrix of the i-th
             * covariate on all edges.
             */
        }
    };

    /* Keep this as is. This will be usefull when symmetric SBM will be
     * implemented */
    struct is_sbm_symmetric
    {
        enum { value=false };
    };

    // parameters
    unsigned int n_parameters;

    /* Here you must put the MODEL model parameters */

    MODEL(SBM & membership, MODEL::network & net)
    {
        /* Here you must intialize the number of parameters (n_parameters) and
         * your model parameters, knowing the membership and the network in the
         * SBM case */
    }
    
    MODEL(LBM & membership, MODEL::network & net)
    {
        /* Here you must intialize the number of parameters (n_parameters) and
         * your model parameters, knowing the membership and the network in the
         * LBM case */
    }

    inline
    Rcpp::List export_to_R()
    {
        Rcpp::List values;
        values["n_parameters"] = n_parameters;

        /* Here you must define the export way of your model parameters in the
         * list values. This is this list which is returned to R. */

        return values;
    }

    inline
    Rcpp::List prediction(SBM & membership, MODEL::network & net)
    {
        Rcpp::List Lpred;

        /* Here you should provide a method to compute the prediction of the
         * network knowing parameters and the membership in the case of SBM. The
         * network object is provided, but only covariate fields shoud be used
         * */

        /* The result is a list for R, as the same organization than the
         * provided network by R. See the constructor of MODEL::network you have
         * previously defined.
         */
        return Lpred;
    }
    
    inline
    Rcpp::List prediction(LBM & membership)
    {
        Rcpp::List Lpred;
        
        /* Here you should provide a method to compute the prediction of the
         * network knowing parameters and the membership in the case of SBM. The
         * network object is provided, but only covariate fields shoud be used
         * */

        /* The result is a list for R, as the same organization than the
         * provided network by R. See the constructor of MODEL::network you have
         * previously defined.
         */
        
        return Lpred;
    }
};

/* Function for the model.
 */

/******************************************************************************
 *
 * Notations :
 *
 * i a row node
 * j a col node
 * q a row class
 * l a col class
 *
 * (for SBM there are no distinctions between row and col nodes and classes)
 *
 * f(i,j,q,l) the prob or density of the variable
 * X_{ij} conditionnally to i belong to class q and j to class l
 *
 * logf(i,j,q,l) is a quicker way to write log(f(i,j,q,l))
 *
 ******************************************************************************/





/* If you have direct formula to compute the following form, use it, else skip
 * this definition and delete this definition.
 *
 * forall i,q ; lZ_{iq} += \sum_{jl} Z_{jl} logf(i,j,q,l)
 *
 * This is usefull for SBM. See bernoulli and poisson for example.
 */

// template<>
// inline
// void e_fixed_step(SBM & membership,
//                   MODEL & model,
//                   MODEL::network & net,
//                   mat & lZ)
// {
// }





/* If you have direct formula to compute the following form, use it, else skip
 * this definition and delete this definition.
 *
 * forall i,q ; lZ1_{iq} += \sum_{jl} Z2_{jl} logf(i,j,q,l)
 * forall j,l ; lZ2_{jl} += \sum_{iq} Z1_{iq} logf(i,j,q,l)
 *
 * This is usefull for LBM. See bernoulli and poisson for example.
 */

// template<>
// inline
// void e_fixed_step(LBM & membership,
//                   MODEL & model,
//                   MODEL::network & net,
//                   mat & lZ1,
//                   mat & lZ2)
// {
// }





/* If you have a explicit formula to compute the maximum of PL in respect to
 * model parameters, knowing the membership, use it, else skip and delete this
 * definition. In the case of SBM.
 *
 * This function must return PL in the maximum. You can use the function PL if
 * you have define one. If you have not define one. A generic function is used,
 * but it is unprobable you have a direct way to compute the maximum of PL
 * without having a direct way to compute PL.
 */

// template<>
// inline
// double m_step(SBM & membership, poisson & model, poisson::network & net)
// {
// }





template<>
inline
double m_step(LBM & membership, poisson & model, poisson::network & net)
{
    model.lambda = (membership.Z1.t() * net.adj * membership.Z2)
                    /
                   (membership.Z1.t() * net.Mones * membership.Z2);

    return
        (
            accu(
                
                    (
                        model.lambda 
                        % 
                        (membership.Z1.t() * net.Mones * membership.Z2)
                    )
                    +
                    (
                        log(model.lambda)
                        %
                        (membership.Z1.t() * net.adj * membership.Z2)
                    )
                )
            + net.accu_log_fact_X
        );
}


