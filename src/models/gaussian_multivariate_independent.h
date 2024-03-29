
class gaussian_multivariate_independent
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
        cube adj;


        /* Here you should add all precomputed values which depends only on the
         * network usefull in various functions, to avoid computing these value
         * many time.
         */
        cube adjZ;
        mat Mones;
        mat MonesZ;
        vec accu_adj_square;
        vec accu_adjZ_square;

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
            Rcpp::List adj_list = network_from_R["adjacency"];
            mat first_mat = Rcpp::as<mat>( adj_list[0] );
            adj.set_size(first_mat.n_rows, first_mat.n_cols, adj_list.size());
            for(int k=0; k<adj_list.size(); k++)
                adj.slice(k) = Rcpp::as<mat>(adj_list[k]);

            Mones = ones<mat>(first_mat.n_rows, first_mat.n_cols);
            MonesZ = fill_diag(Mones,0);
            adjZ.set_size(adj.n_rows,adj.n_cols,adj.n_slices);
            accu_adj_square.set_size(adj.n_slices);
            accu_adjZ_square.set_size(adj.n_slices);
            for(unsigned int k=0; k<adj.n_slices; k++)
            {
                adjZ.slice(k) = fill_diag(adj.slice(k),0);
                accu_adj_square(k) = accu(adj.slice(k) % adj.slice(k));
                accu_adjZ_square(k) = accu(adjZ.slice(k) % adjZ.slice(k));
            }
        }
    };

    // parameters
    unsigned int n_parameters;

    /* Here you must put the gaussian_multivariate_independent model parameters */
    cube mu;
    vec sigma2;

    gaussian_multivariate_independent(SBM & membership, gaussian_multivariate_independent::network & net)
    {
        /* Here you must intialize the number of parameters (n_parameters) and
         * your model parameters, knowing the membership and the network in the
         * SBM case
         */
        n_parameters = membership.Z.n_cols * membership.Z.n_cols * net.adj.n_slices + net.adj.n_slices;
        mu.set_size(membership.Z.n_cols, membership.Z.n_cols, net.adj.n_slices);
        sigma2.set_size(net.adj.n_slices);
    }
    
    gaussian_multivariate_independent(SBM_sym & membership, gaussian_multivariate_independent::network & net)
    {
        /* Here you must intialize the number of parameters (n_parameters) and
         * your model parameters, knowing the membership and the network in the
         * SBM_sym case
         */
        n_parameters = membership.Z.n_cols * (membership.Z.n_cols+1)/2 * net.adj.n_slices + net.adj.n_slices;
        mu.set_size(membership.Z.n_cols, membership.Z.n_cols, net.adj.n_slices);
        sigma2.set_size(net.adj.n_slices);
    }
    
    gaussian_multivariate_independent(LBM & membership, gaussian_multivariate_independent::network & net)
    {
        /* Here you must intialize the number of parameters (n_parameters) and
         * your model parameters, knowing the membership and the network in the
         * LBM case 
         */
        n_parameters = membership.Z1.n_cols * membership.Z2.n_cols * net.adj.n_slices + net.adj.n_slices;
        mu.set_size(membership.Z1.n_cols, membership.Z2.n_cols, net.adj.n_slices);
        sigma2.set_size(net.adj.n_slices);
    }

    inline
    Rcpp::List export_to_R()
    {
        Rcpp::List values;
        values["n_parameters"] = n_parameters;
     
        values["sigma2"] = sigma2;
        Rcpp::List mu_as_list(mu.n_slices);
        for(unsigned int k=0;k<mu.n_slices;k++)
        {
            mu_as_list[k]=mu.slice(k);
        }
        values["mu"] = mu_as_list;

        /* Here you must define the export way of your model parameters in the
         * list values. This is this list which is returned to R. */

        return values;
    }

    /* Keep this. When this is not usefull, this is not used. */
    inline vec to_vector();
    gaussian_multivariate_independent(SBM &, const vec &);
    gaussian_multivariate_independent(SBM_sym &, const vec &);
    gaussian_multivariate_independent(LBM &, const vec &);
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





/* Usefull, Optional
 *
 * If you have direct formula to compute the following form, use it, else skip
 * this definition and delete this definition.
 *
 * forall i,q ; lZ_{iq} += \sum_{jl} Z_{jl} (logf(i,j,q,l) + logf(j,i,l,q)
 *
 * This is usefull for SBM. See bernoulli and poisson for example.
 */

template<>
inline
void e_fixed_step(SBM & membership,
                  gaussian_multivariate_independent & model,
                  gaussian_multivariate_independent::network & net,
                  mat & lZ)
{
    for(unsigned int k=0;k<net.adj.n_slices;k++)
    {
        lZ += 1.0/(2*model.sigma2(k)) * (
                - net.MonesZ * membership.Z * (model.mu.slice(k).t() % model.mu.slice(k).t())
                + 2 * net.adjZ.slice(k) * membership.Z * model.mu.slice(k).t()
                - net.MonesZ.t() * membership.Z * (model.mu.slice(k) % model.mu.slice(k))
                + 2 * net.adjZ.slice(k).t() * membership.Z * model.mu.slice(k));
    }
}





/* Usefull, Optional
 *
 * If you have direct formula to compute the following form, use it, else skip
 * this definition and delete this definition.
 *
 * forall i,q ; lZ_{iq} += \sum_{jl} Z_{jl} logf(i,j,q,l)
 *
 * This is usefull for SBM_sym. See bernoulli and poisson for example.
 */

template<>
inline
void e_fixed_step(SBM_sym & membership,
                  gaussian_multivariate_independent & model,
                  gaussian_multivariate_independent::network & net,
                  mat & lZ)
{
    for(unsigned int k=0;k<net.adj.n_slices;k++)
    {
        lZ += 1.0/(2*model.sigma2(k)) * (
                - net.MonesZ * membership.Z * (model.mu.slice(k).t() % model.mu.slice(k).t())
                + 2 * net.adjZ.slice(k) * membership.Z * model.mu.slice(k).t());
    }
}





/* Usefull, Optional
 *
 * If you have direct formula to compute the following form, use it, else skip
 * this definition and delete this definition.
 *
 * forall i,q ; lZ1_{iq} += \sum_{jl} Z2_{jl} logf(i,j,q,l)
 * forall j,l ; lZ2_{jl} += \sum_{iq} Z1_{iq} logf(i,j,q,l)
 *
 * This is usefull for LBM. See bernoulli and poisson for example.
 */

template<>
inline
void e_fixed_step(LBM & membership,
                  gaussian_multivariate_independent & model,
                  gaussian_multivariate_independent::network & net,
                  mat & lZ1,
                  mat & lZ2)
{
    for(unsigned int k=0;k<net.adj.n_slices;k++)
    {
        lZ1 += 1.0/(2*model.sigma2(k)) * (
                - net.Mones * membership.Z2 * (model.mu.slice(k).t() % model.mu.slice(k).t())
                + 2 * net.adj.slice(k) * membership.Z2 * model.mu.slice(k).t());
        lZ2 += 1.0/(2*model.sigma2(k)) * (
                - net.Mones.t() * membership.Z1 * (model.mu.slice(k) % model.mu.slice(k))
                + 2 * net.adj.slice(k).t() * membership.Z1 * model.mu.slice(k));
    }
}





/* Usefull, Optional
 *
 * If you have a explicit formula to compute the maximum of PL in respect to
 * model parameters, knowing the membership, use it, else skip and delete this
 * definition. In the case of SBM.
 *
 * This function must return PL in the maximum.
 * 
 * This is usefull for SBM. See bernoulli and poisson for example.
 */

template<>
inline
double m_step(SBM & membership,
              gaussian_multivariate_independent & model,
              gaussian_multivariate_independent::network & net)
{
    mat provdiv = membership.Z.t() * net.MonesZ * membership.Z;
    for(unsigned int k=0; k<net.adj.n_slices;k++)
    {
        model.mu.slice(k) = (membership.Z.t() * net.adjZ.slice(k) * membership.Z)
                            /
                            provdiv;
    }

    vec all_accu = zeros<vec>(net.adj.n_slices);
    for(unsigned int k=0;k<net.adj.n_slices;k++)
    {
        all_accu(k) = net.accu_adjZ_square(k) +
            accu(
                (
                    (model.mu.slice(k) % model.mu.slice(k))
                    %
                    ( membership.Z.t() * net.MonesZ * membership.Z )
                )
                -
                (
                    2 * model.mu.slice(k)
                    %
                    (membership.Z.t() * net.adjZ.slice(k) * membership.Z)
                )
            );

        model.sigma2(k) = 1.0/(membership.Z.n_rows * (membership.Z.n_rows-1)) * all_accu(k);
    }


    double PL= -.5*(membership.Z.n_rows * (membership.Z.n_rows-1))*accu(log(2*M_PI*model.sigma2));
    for(unsigned int k=0;k<net.adj.n_slices;k++)
        PL+= -1.0/(2*model.sigma2(k))*all_accu(k);

    return PL;
}





/* Usefull, Optional
 *
 * If you have a explicit formula to compute the maximum of PL in respect to
 * model parameters, knowing the membership, use it, else skip and delete this
 * definition. In the case of SBM_sym.
 *
 * This function must return PL in the maximum.
 * 
 * This is usefull for SBM_sym. See bernoulli and poisson for example.
 */

template<>
inline
double m_step(SBM_sym & membership,
              gaussian_multivariate_independent & model,
              gaussian_multivariate_independent::network & net)
{
    return m_step<SBM>(membership,model,net)/2;
}





/* Usefull, Optional
 *
 * If you have a explicit formula to compute the maximum of PL in respect to
 * model parameters, knowing the membership, use it, else skip and delete this
 * definition. In the case of LBM.
 *
 * This function must return PL in the maximum.
 * 
 * This is usefull for LBM. See bernoulli and poisson for example.
 */

template<>
inline
double m_step(LBM & membership,
              gaussian_multivariate_independent & model,
              gaussian_multivariate_independent::network & net)
{
    mat provdiv = membership.Z1.t() * net.Mones * membership.Z2;
    for(unsigned int k=0; k<net.adj.n_slices;k++)
    {
        model.mu.slice(k) = (membership.Z1.t() * net.adj.slice(k) * membership.Z2)
                            /
                            provdiv;
    }

    vec all_accu = zeros<vec>(net.adj.n_slices);
    for(unsigned int k=0;k<net.adj.n_slices;k++)
    {
        all_accu(k) = net.accu_adj_square(k) +
            accu(
                (
                    (model.mu.slice(k) % model.mu.slice(k))
                    %
                    ( membership.Z1.t() * net.Mones * membership.Z2 )
                )
                -
                (
                    2 * model.mu.slice(k)
                    %
                    (membership.Z1.t() * net.adj.slice(k) * membership.Z2)
                )
            );
            model.sigma2(k) = 1.0/(membership.Z1.n_rows * membership.Z2.n_rows) * all_accu(k);
    }


    double PL =  -.5*(membership.Z1.n_rows * membership.Z2.n_rows)*accu(log(2*M_PI*model.sigma2));
    for(unsigned int k=0;k<net.adj.n_slices;k++) 
        PL +=  -1.0/(2*model.sigma2(k))*all_accu(k);

    return PL;
}





/* If you have defined m_step(SBM &, gaussian_multivariate_independent &, gaussian_multivariate_independent::network &)
 * Then
 *     Useless
 * Else
 *     Usefull, Optional
 *
 * If you have direct formula to compute PL, without loops over (i,j,q,l), in
 * the case of SBM, use it.
 *
 * See poisson_covariate for example.
 */

// template<>
// inline
// double PL(gaussian_multivariate_independent & model,
//           SBM & membership,
//           gaussian_multivariate_independent::network & net)





/* If you have defined m_step(SBM_sym &, gaussian_multivariate_independent &, gaussian_multivariate_independent::network &)
 * Then
 *     Useless
 * Else
 *     Usefull, Optional
 *
 * If you have direct formula to compute PL, without loops over (i,j,q,l), in
 * the case of SBM_sym, use it.
 *
 * See poisson_covariate for example.
 */

// template<>
// inline
// double PL(gaussian_multivariate_independent & model,
//           SBM_sym & membership,
//           gaussian_multivariate_independent::network & net)





/* If you have defined m_step(LBM &, gaussian_multivariate_independent &, gaussian_multivariate_independent::network &)
 * Then
 *     Useless
 * Else
 *     Usefull, Optional
 *
 * If you have direct formula to compute PL, without loops over (i,j,q,l), in
 * the case of LBM, use it.
 *
 * See poisson_covariate for example.
 */

// template<>
// inline
// double PL(gaussian_multivariate_independent & model,
//           LBM & membership,
//           gaussian_multivariate_independent::network & net)





/* If one of the following specialization
 *     - m_step(LBM &, gaussian_multivariate_independent &, gaussian_multivariate_independent::network)
 *     - m_step(SBM &, gaussian_multivariate_independent &, gaussian_multivariate_independent::network)
 *     - m_step(SBM_sym &, gaussian_multivariate_independent &, gaussian_multivariate_independent::network)
 *   is not defined
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 *
 * Export the model parameters as a vector.
 *
 * See naive_bernoulli and poisson_covariates for example.
 */

// inline
// vec gaussian_multivariate_independent::to_vector()
// {
// }





/* If m_step(SBM &, gaussian_multivariate_independent &, gaussian_multivariate_independent::network) is not defined
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 * 
 * Constructor of the model, from a vector of parameter named vectorized with
 * the given membership of type SBM.
 *
 * See naive_bernoulli for example.
 */

// gaussian_multivariate_independent::gaussian_multivariate_independent(SBM & membership, const vec & vectorized)
// {
// }





/* If m_step(SBM_sym &, gaussian_multivariate_independent &, gaussian_multivariate_independent::network) is not defined
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 * 
 * Constructor of the model, from a vector of parameter named vectorized with
 * the given membership of type SBM_sym.
 *
 * See naive_bernoulli for example.
 */

// gaussian_multivariate_independent::gaussian_multivariate_independent(SBM_sym & membership, const vec & vectorized)
// {
// }





/* If m_step(LBM &, gaussian_multivariate_independent &, gaussian_multivariate_independent::network) is not defined
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 * 
 * Constructor of the model, from a vector of parameter named vectorized with
 * the given membership of type LBM.
 *
 * See naive_bernoulli for example.
 */

// gaussian_multivariate_independent::gaussian_multivariate_independent(LBM & membership, const vec & vectorized)
// {
// }





/* If one of the following specialization
 *     - m_step(LBM &, gaussian_multivariate_independent &, gaussian_multivariate_independent::network)
 *     - m_step(SBM &, gaussian_multivariate_independent &, gaussian_multivariate_independent::network)
 *     - m_step(SBM_sym &, gaussian_multivariate_independent &, gaussian_multivariate_independent::network)
 *   is not defined
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 *
 * From the given model, in the given direction for model parameters, return the
 * maximum betwenn 1 and the maximum step the descend algorithm can do.
 *
 * See naive_bernoulli for example.
 *
 * If you don't have any constraint on you parameter vector, return 1, else, use
 * your brain.
 *
 * You can specialize this template in three function for SBM, LBM, and SBM_sym
 * if you need
 */

// template<class membership_type>
// inline
// double maximum_step_in_direction(membership_type & membership,
//                                  gaussian_multivariate_independent & model,
//                                  gaussian_multivariate_independent::network & net,
//                                  vec & direction)
// {
// }





/* If you have defined m_step(SBM &, gaussian_multivariate_independent &, gaussian_multivariate_independent::network &)
 * Then
 *     Useless
 * Else
 *     Usefull, Optional
 *
 * If you have a direct formula to compute the gradient of PL in respect to your
 * parameters vectors without loops over (i,j,q,l), in the case of SBM, use it.
 *
 * See poisson_covariates for example.
 */

// template<>
// inline
// vec grad(gaussian_multivariate_independent & model,
//          SBM & membership,
//          gaussian_multivariate_independent::network & net)
// {
// }





/* If you have defined m_step(SBM_sym &, gaussian_multivariate_independent &, gaussian_multivariate_independent::network &)
 * Then
 *     Useless
 * Else
 *     Usefull, Optional
 *
 * If you have a direct formula to compute the gradient of PL in respect to your
 * parameters vectors without loops over (i,j,q,l), in the case of SBM, use it.
 *
 * See poisson_covariates for example.
 */

// template<>
// inline
// vec grad(gaussian_multivariate_independent & model,
//          SBM_sym & membership,
//          gaussian_multivariate_independent::network & net)
// {
// }





/* If you have defined m_step(LBM &, gaussian_multivariate_independent &, gaussian_multivariate_independent::network &)
 * Then
 *     Useless
 * Else
 *     Usefull, Optional
 *
 * If you have a direct formula to compute the gradient of PL in respect to your
 * parameters vectors without loops over (i,j,q,l), in the case of LBM, use it.
 *
 * See poisson_covariates for example.
 */

// template<>
// inline
// vec grad(gaussian_multivariate_independent & model,
//          LBM & membership,
//          gaussian_multivariate_independent::network & net)
// {
// }





/* If one of the following specialization
 *     - grad(gaussian_multivariate_independent &, LBM &, gaussian_multivariate_independent::network)
 *     - grad(gaussian_multivariate_independent &, SBM &, gaussian_multivariate_independent::network)
 *     - grad(gaussian_multivariate_independent &, SBM_sym &, gaussian_multivariate_independent::network)
 *   is Usefull and not defined
 * Then
 *     Usefull, Optional
 * Else
 *     Useless
 *
 * The gradient of logf(i,j,q,l) in respect to the parameter vector.
 */

// template<>
// inline
// vec grad_logf(gaussian_multivariate_independent & model,
//               gaussian_multivariate_independent::network & net,
//               unsigned int i,
//               unsigned int j,
//               unsigned int q,
//               unsigned int l)
// {
// }





/* If grad_logf(gaussian_multivariate_independent &,
 *         gaussian_multivariate_independent::network &,
 *         unsigned int,
 *         unsigned int,
 *         unsigned int,
 *         unsigned int)
 *    is Usefull and not defined
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 *
 *
 * The derivative of logf(i,j,q,l) in respect to the k-th parameters of the
 * parameter vector.
 *
 * See naive_bernoulli for example.
 */

// inline
// double grad_logf(gaussian_multivariate_independent & model,
//                  gaussian_multivariate_independent::network & net,
//                  unsigned int i,
//                  unsigned int j,
//                  unsigned int q,
//                  unsigned int l,
//                  unsigned int k)
// {
// }





/* If one of the following specialization
 *     - PL(gaussian_multivariate_independent &, SBM &, gaussian_multivariate_independent::network)
 *     - PL(gaussian_multivariate_independent &, SBM_sym &, gaussian_multivariate_independent::network)
 *     - PL(gaussian_multivariate_independent &, LBM &, gaussian_multivariate_independent::network)
 *   is Usefull and not definied
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 *
 * logf(i,j,q,l) see notations.
 */

// inline
// double logf(gaussian_multivariate_independent & model,
//             gaussian_multivariate_independent::network & net,
//             unsigned int i,
//             unsigned int j,
//             unsigned int q,
//             unsigned int l)
// {
// }
