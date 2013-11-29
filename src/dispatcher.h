
#include "conf.h"

#include <RcppArmadillo.h>
#include <string>

using namespace arma;

#include "misc_functions.h"
#include "memberships/sbm.h"
#include "memberships/lbm.h"
#include "generic_functions.h"
#include "model_base.h"
#include "general_em.h"
#include "models/naive_bernoulli.h"
#include "models/bernoulli.h"


template<class membership_type, class model_type, class network_type>
Rcpp::List estim(membership_type & membership_init,
           Rcpp::List & network_from_R)
{
    network_type net(network_from_R);

    result<membership_type, model_type> res =
        em<membership_type, model_type>( membership_init, net );

    return res.export_to_R();
}

template<class membership_type>
Rcpp::List dispatcher_model(membership_type & membership_init,
                      std::string & model_name,
                      Rcpp::List & network_from_R)
{

    /* Each model must have a entry here */

    if(model_name == "naive_bernoulli")
        return estim<membership_type, naive_bernoulli, naive_bernoulli::network>(
                membership_init, network_from_R);

    if(model_name == "bernoulli")
    {

        return estim<membership_type, bernoulli, bernoulli::network>(
                membership_init, network_from_R);
    }

    // What we are doing here ? Nothing.
    return Rcpp::List();
}

template<class membership_type>
Rcpp::List init_membership_and_dispatcher_model(Rcpp::List & membership_from_R,
                                          std::string & model_name,
                                          Rcpp::List & network_from_R)
{
    membership_type membership_init(membership_from_R);
    return dispatcher_model(membership_init, model_name, network_from_R);
}


Rcpp::List distpatcher_membership_model(std::string & membership_name,
                                  Rcpp::List & membership_init_from_R,
                                  std::string & model_name,
                                  Rcpp::List & network_from_R)
{
    /* Each membership must have a entry here */

    if(membership_name == "SBM")
        return init_membership_and_dispatcher_model<SBM>(membership_init_from_R,
                                                         model_name,
                                                         network_from_R);
    
    if(membership_name == "LBM")
        return init_membership_and_dispatcher_model<LBM>(membership_init_from_R,
                                                         model_name,
                                                         network_from_R);

    // What we are doing here ? Nothing.
    return Rcpp::List();

}
                                          

RcppExport
SEXP dispatcher(SEXP membership_name,
                SEXP membership_init_from_R,
                SEXP model_name,
                SEXP network_from_R);

