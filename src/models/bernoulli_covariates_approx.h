
/*
 * This is package skeleton for bernoulli_covariates_approx model
 *
 * You should delete this header when the bernoulli_covariates_approx model is ready for use to avoid
 * any confusion
 */

class bernoulli_covariates_approx
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

        mat adj;
        cube covariates;

        mat adjZD;
        mat Mones;
        mat Monest;
        mat MonesZD;
        mat adjt;


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
            adj = Rcpp::as<mat>(network_from_R["adjacency"]);

            Rcpp::List covariates_list = network_from_R["covariates"];
            covariates.set_size(adj.n_rows,adj.n_cols,covariates_list.size());
            for(int k=0; k<covariates_list.size(); k++)
                covariates.slice(k) = Rcpp::as<mat>(covariates_list[k]);

            adjZD = fill_diag(adj,0);
            Mones = ones<mat>(adj.n_rows,adj.n_cols);
            MonesZD = fill_diag(Mones,0);
            adjt = adj.t();
            Monest = Mones.t();

        }
    };

    // parameters
    unsigned int n_parameters;
    bool symmetric;

    /* Here you must put the bernoulli_covariates_approx model parameters */
    mat m;
    colvec beta;

    bernoulli_covariates_approx(SBM & membership, bernoulli_covariates_approx::network & net)
    {
        /* Here you must intialize the number of parameters (n_parameters) and
         * your model parameters, knowing the membership and the network in the
         * SBM case
         */
        m = (membership.Z.t() * net.adjZD * membership.Z)
            /
            (membership.Z.t() * net.MonesZD * membership.Z);

        m=logit(m);

        beta.set_size(net.covariates.n_slices);
        beta.fill(0);

        n_parameters = m.n_elem + beta.n_elem;
        symmetric=false;
    }
    
    bernoulli_covariates_approx(SBM_sym & membership, bernoulli_covariates_approx::network & net)
    {
        /* Here you must intialize the number of parameters (n_parameters) and
         * your model parameters, knowing the membership and the network in the
         * SBM_sym case
         */
        m = (membership.Z.t() * net.adjZD * membership.Z)
            /
            (membership.Z.t() * net.MonesZD * membership.Z);

        m=logit(m);

        beta.set_size(net.covariates.n_slices);
        beta.fill(0);

        n_parameters = m.n_rows*(m.n_rows+1)/2 + beta.n_elem;
        symmetric=true;
    }
    
    bernoulli_covariates_approx(LBM & membership, bernoulli_covariates_approx::network & net)
    {
        /* Here you must intialize the number of parameters (n_parameters) and
         * your model parameters, knowing the membership and the network in the
         * LBM case 
         */
        m = (membership.Z1.t() * net.adj * membership.Z2)
            /
            (membership.Z1.t() * net.Mones * membership.Z2);

        m=logit(m);

        beta.set_size(net.covariates.n_slices);
        beta.fill(0);

        n_parameters = m.n_rows * m.n_cols + net.covariates.n_slices;
        symmetric=false;
    }

    inline
    Rcpp::List export_to_R()
    {
        Rcpp::List values;
        values["n_parameters"] = n_parameters;

        /* Here you must define the export way of your model parameters in the
         * list values. This is this list which is returned to R. */

        values["m"] = m;
        values["beta"] = beta;

        return values;
    }

    /* Keep this. When this is not usefull, this is not used. */
    inline vec to_vector();
    bernoulli_covariates_approx(SBM &, const vec &);
    bernoulli_covariates_approx(SBM_sym &, const vec &);
    bernoulli_covariates_approx(LBM &, const vec &);
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
                  bernoulli_covariates_approx & model,
                  bernoulli_covariates_approx::network & net,
                  mat & lZ)
{
    mat B = fill_diag(compute_B(model.beta, net.covariates),0);
    mat BhpB = B%B;
    mat mm = -model.m;
    mat smm = sigmo(mm);
    mat lsmm = log(smm);
    mat sm = sigmo(model.m);
    mat smhpsmm = sm % smm;

    lZ += net.adjZD * membership.Z * model.m.t()
        + net.MonesZD * membership.Z * lsmm.t()
        - B * membership.Z * sm.t()
        - .5 * BhpB * membership.Z * smhpsmm.t()
        + net.adjZD.t() * membership.Z * model.m
        + net.MonesZD * membership.Z * lsmm
        - B.t() * membership.Z * sm
        - .5 * BhpB.t() * membership.Z * smhpsmm;
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
                  bernoulli_covariates_approx & model,
                  bernoulli_covariates_approx::network & net,
                  mat & lZ)
{
    mat B = fill_diag(compute_B(model.beta, net.covariates),0);
    mat BhpB = B%B;
    mat mm = -model.m;
    mat smm = sigmo(mm);
    mat lsmm = log(smm);
    mat sm = sigmo(model.m);
    mat smhpsmm = sm % smm;

    lZ += net.adjZD * membership.Z * model.m
        + net.MonesZD * membership.Z * lsmm
        - B * membership.Z * sm
        - .5 * BhpB * membership.Z * smhpsmm;
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
                  bernoulli_covariates_approx & model,
                  bernoulli_covariates_approx::network & net,
                  mat & lZ1,
                  mat & lZ2)
{
    mat B = compute_B(model.beta, net.covariates);
    mat BhpB = B%B;
    mat mm = -model.m;
    mat smm = sigmo(mm);
    mat lsmm = log(smm);
    mat sm = sigmo(model.m);
    mat smhpsmm = sm % smm;

    lZ1 += net.adj * membership.Z2 * model.m.t()
        + net.Mones * membership.Z2 * lsmm.t()
        - B * membership.Z2 * sm.t()
        - .5 * BhpB * membership.Z2 * smhpsmm.t();
    
    lZ2 += net.adjt * membership.Z1 * model.m
        + net.Monest * membership.Z1 * lsmm
        - B.t() * membership.Z1 * sm
        - .5 * BhpB.t() * membership.Z1 * smhpsmm;
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

// template<>
// inline
// double m_step(SBM & membership,
//               bernoulli_covariates_approx & model,
//               bernoulli_covariates_approx::network & net)
// {
// }





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

// template<>
// inline
// double m_step(SBM_sym & membership,
//               bernoulli_covariates_approx & model,
//               bernoulli_covariates_approx::network & net)
// {
// }





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

// template<>
// inline
// double m_step(LBM & membership,
//               bernoulli_covariates_approx & model,
//               bernoulli_covariates_approx::network & net)
// {
// }





/* If you have defined m_step(SBM &, bernoulli_covariates_approx &, bernoulli_covariates_approx::network &)
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

template<>
inline
double PL(bernoulli_covariates_approx & model,
          SBM & membership,
          bernoulli_covariates_approx::network & net)
{
    mat B = fill_diag(compute_B(model.beta, net.covariates),0);
    mat mm = -model.m;
    return(
            accu(
                (membership.Z.t() * net.adjZD * membership.Z) % model.m
                )
            +
            accu(
                net.adjZD % B
                )
            +
            accu(
                membership.Z * log(sigmo(mm)) * membership.Z.t()
                )
            -
            accu(
                ( membership.Z.t() * B * membership.Z )
                %
                sigmo(model.m)
                )
            - .5 *
            accu(
                ( membership.Z.t() * (B%B) * membership.Z )
                %
                sigmo(model.m) % sigmo(mm)
                )
          );
}





/* If you have defined m_step(SBM_sym &, bernoulli_covariates_approx &, bernoulli_covariates_approx::network &)
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

template<>
inline
double PL(bernoulli_covariates_approx & model,
          SBM_sym & membership,
          bernoulli_covariates_approx::network & net)
{
    mat B = fill_diag(compute_B(model.beta, net.covariates),0);
    mat mm = -model.m;
    return( .5*(
            accu(
                (membership.Z.t() * net.adjZD * membership.Z) % model.m
                )
            +
            accu(
                net.adjZD % B
                )
            +
            accu(
                membership.Z * log(sigmo(mm)) * membership.Z.t()
                )
            -
            accu(
                ( membership.Z.t() * B * membership.Z )
                %
                sigmo(model.m)
                )
            - .5 *
            accu(
                ( membership.Z.t() * (B%B) * membership.Z )
                %
                sigmo(model.m) % sigmo(mm)
                )
          ) );
}





/* If you have defined m_step(LBM &, bernoulli_covariates_approx &, bernoulli_covariates_approx::network &)
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

template<>
inline
double PL(bernoulli_covariates_approx & model,
          LBM & membership,
          bernoulli_covariates_approx::network & net)
{
    mat B = compute_B(model.beta, net.covariates);
    mat mm = -model.m;
    return(
            accu(
                (membership.Z1.t() * net.adj * membership.Z2) % model.m
                )
            +
            accu(
                net.adj % B
                )
            +
            accu(
                membership.Z1 * log(sigmo(mm)) * membership.Z2.t()
                )
            -
            accu(
                ( membership.Z1.t() * B * membership.Z2 )
                %
                sigmo(model.m)
                )
            - .5 *
            accu(
                ( membership.Z1.t() * (B%B) * membership.Z2 )
                %
                sigmo(model.m) % sigmo(mm)
                )
          );
}





/* If one of the following specialization
 *     - m_step(LBM &, bernoulli_covariates_approx &, bernoulli_covariates_approx::network)
 *     - m_step(SBM &, bernoulli_covariates_approx &, bernoulli_covariates_approx::network)
 *     - m_step(SBM_sym &, bernoulli_covariates_approx &, bernoulli_covariates_approx::network)
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

inline
vec bernoulli_covariates_approx::to_vector()
{
    vec out(n_parameters);
    vec vm = (symmetric) ? vech(m) : reshape(m,m.n_elem,1);

    out.subvec(0,vm.n_elem-1) = vm;
    out.subvec(vm.n_elem,n_parameters-1) = beta;
    return(out);
}





/* If m_step(SBM &, bernoulli_covariates_approx &, bernoulli_covariates_approx::network) is not defined
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

bernoulli_covariates_approx::bernoulli_covariates_approx(SBM & membership, const vec & vectorized)
{
    unsigned int Q=membership.Z.n_cols;
    m = reshape(vectorized.subvec(0,Q*Q-1),Q,Q);
    beta = vectorized.subvec(Q*Q,vectorized.n_elem-1);

    n_parameters = vectorized.n_elem;
    symmetric = false;
}





/* If m_step(SBM_sym &, bernoulli_covariates_approx &, bernoulli_covariates_approx::network) is not defined
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

bernoulli_covariates_approx::bernoulli_covariates_approx(SBM_sym & membership, const vec & vectorized)
{
    unsigned int Q=membership.Z.n_cols;
    m = unvech(vectorized.subvec(0,Q*(Q+1)/2-1));
    beta = vectorized.subvec(Q*(Q+1)/2,vectorized.n_elem-1);

    n_parameters = vectorized.n_elem;
    symmetric = true;
}





/* If m_step(LBM &, bernoulli_covariates_approx &, bernoulli_covariates_approx::network) is not defined
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

bernoulli_covariates_approx::bernoulli_covariates_approx(LBM & membership, const vec & vectorized)
{
    unsigned int Q1=membership.Z1.n_cols;
    unsigned int Q2=membership.Z2.n_cols;
    m = reshape(vectorized.subvec(0,Q1*Q2-1),Q1,Q2);
    beta = vectorized.subvec(Q1*Q2,vectorized.n_elem-1);

    n_parameters = vectorized.n_elem;
    symmetric = false;
}





/* If one of the following specialization
 *     - m_step(LBM &, bernoulli_covariates_approx &, bernoulli_covariates_approx::network)
 *     - m_step(SBM &, bernoulli_covariates_approx &, bernoulli_covariates_approx::network)
 *     - m_step(SBM_sym &, bernoulli_covariates_approx &, bernoulli_covariates_approx::network)
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
 */

inline
double maximum_step_in_direction(bernoulli_covariates_approx & model, vec & direction)
{
    return(1);
}





/* If you have defined m_step(SBM &, bernoulli_covariates_approx &, bernoulli_covariates_approx::network &)
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

template<>
inline
vec grad(bernoulli_covariates_approx & model,
         SBM & membership,
         bernoulli_covariates_approx::network & net)
{
    mat B = fill_diag(compute_B(model.beta,net.covariates),0);
    mat mm = -model.m;
    mat sm = sigmo(model.m);
    mat smm = sigmo(mm);

    mat dm = membership.Z.t() * net.adjZD * membership.Z
           - (membership.Z.t() * net.MonesZD * membership.Z) % sm
           + (membership.Z.t() * B * membership.Z) % sm % smm
           - .5 * (membership.Z.t() * (B % B) * membership.Z) % sm % smm % (smm - sm);

    mat M = net.adjZD 
          - fill_diag(
                  membership.Z * sm * membership.Z.t()
                  ,0)
          - (membership.Z * (sm % smm) * membership.Z.t()) % B;

    vec dbeta(model.beta.n_elem);
    for(unsigned int k=0;k<dbeta.n_elem;k++)
    {
        dbeta(k) = accu( M % net.covariates.slice(k));
    }

    vec out(model.n_parameters);
    vec vdm = reshape(dm,dm.n_elem,1);

    out.subvec(0,vdm.n_elem-1) = vdm;
    out.subvec(vdm.n_elem,model.n_parameters-1) = dbeta;

    return(out);
}





/* If you have defined m_step(SBM_sym &, bernoulli_covariates_approx &, bernoulli_covariates_approx::network &)
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

template<>
inline
vec grad(bernoulli_covariates_approx & model,
         SBM_sym & membership,
         bernoulli_covariates_approx::network & net)
{
    mat B = fill_diag(compute_B(model.beta,net.covariates),0);
    mat mm = -model.m;
    mat sm = sigmo(model.m);
    mat smm = sigmo(mm);

    mat dm = membership.Z.t() * net.adjZD * membership.Z
           - (membership.Z.t() * net.MonesZD * membership.Z) % sm
           + (membership.Z.t() * B * membership.Z) % sm % smm
           - .5 * (membership.Z.t() * (B % B) * membership.Z) % sm % smm % (smm - sm);

    mat M = net.adjZD 
          - fill_diag(
                  membership.Z * sm * membership.Z.t()
                  ,0)
          - (membership.Z * (sm % smm) * membership.Z.t()) % B;

    vec dbeta(model.beta.n_elem);
    for(unsigned int k=0;k<dbeta.n_elem;k++)
    {
        dbeta(k) = accu( M % net.covariates.slice(k));
    }

    vec out(model.n_parameters);
    vec vdm = vech(dm);

    out.subvec(0,vdm.n_elem-1) = vdm;
    out.subvec(vdm.n_elem,model.n_parameters-1) = dbeta;

    return(.5*out);
}





/* If you have defined m_step(LBM &, bernoulli_covariates_approx &, bernoulli_covariates_approx::network &)
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

template<>
inline
vec grad(bernoulli_covariates_approx & model,
         LBM & membership,
         bernoulli_covariates_approx::network & net)
{
    mat B = compute_B(model.beta,net.covariates);
    mat mm = -model.m;
    mat sm = sigmo(model.m);
    mat smm = sigmo(mm);

    mat dm = membership.Z1.t() * net.adj * membership.Z2
           - (membership.Z1.t() * net.Mones * membership.Z2) % sm
           + (membership.Z1.t() * B * membership.Z2) % sm % smm
           - .5 * (membership.Z1.t() * (B % B) * membership.Z2) % sm % smm % (smm - sm);

    mat M = net.adj
          - membership.Z1 * sm * membership.Z2.t()
          - (membership.Z1 * (sm % smm) * membership.Z2.t()) % B;

    vec dbeta(model.beta.n_elem);
    for(unsigned int k=0;k<dbeta.n_elem;k++)
    {
        dbeta(k) = accu( M % net.covariates.slice(k));
    }

    vec out(model.n_parameters);
    vec vdm = reshape(dm,dm.n_elem,1);

    out.subvec(0,vdm.n_elem-1) = vdm;
    out.subvec(vdm.n_elem,model.n_parameters-1) = dbeta;

    return(out);
}





/* If one of the following specialization
 *     - grad(bernoulli_covariates_approx &, LBM &, bernoulli_covariates_approx::network)
 *     - grad(bernoulli_covariates_approx &, SBM &, bernoulli_covariates_approx::network)
 *     - grad(bernoulli_covariates_approx &, SBM_sym &, bernoulli_covariates_approx::network)
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
// vec grad_logf(bernoulli_covariates_approx & model,
//               bernoulli_covariates_approx::network & net,
//               unsigned int i,
//               unsigned int j,
//               unsigned int q,
//               unsigned int l)
// {
// }





/* If grad_logf(bernoulli_covariates_approx &,
 *         bernoulli_covariates_approx::network &,
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
// double grad_logf(bernoulli_covariates_approx & model,
//                  bernoulli_covariates_approx::network & net,
//                  unsigned int i,
//                  unsigned int j,
//                  unsigned int q,
//                  unsigned int l,
//                  unsigned int k)
// {
// }





/* If one of the following specialization
 *     - PL(bernoulli_covariates_approx &, SBM &, bernoulli_covariates_approx::network)
 *     - PL(bernoulli_covariates_approx &, SBM_sym &, bernoulli_covariates_approx::network)
 *     - PL(bernoulli_covariates_approx &, LBM &, bernoulli_covariates_approx::network)
 *   is Usefull and not definied
 * Then
 *     Usefull, Mandatory
 * Else
 *     Useless
 *
 * logf(i,j,q,l) see notations.
 */

// inline
// double logf(bernoulli_covariates_approx & model,
//             bernoulli_covariates_approx::network & net,
//             unsigned int i,
//             unsigned int j,
//             unsigned int q,
//             unsigned int l)
// {
// }

template<>
inline
void print_model(bernoulli_covariates_approx & model)
{
    model.m.print("m:");
    model.beta.print("beta:");
}

