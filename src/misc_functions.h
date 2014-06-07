
template<class T>
inline
void boundaries(T & obj, double b1, double b2)
{
    for(unsigned int i=0; i<obj.n_elem; i++)
    {
        if(obj(i)<b1)
            obj(i)=b1;
        else
            if(obj(i)>b2)
                obj(i)=b2;
    }
}

inline
mat fill_diag(mat the_mat, double value)
{
    unsigned int n=((the_mat.n_rows<the_mat.n_cols)?the_mat.n_rows:the_mat.n_cols);
    for(unsigned int i=0; i<n; i++)
        the_mat(i,i)=value;
    return the_mat;
}

inline
void accu_log_fact(mat & M, double & with_diag, double & without_diag)
{
    with_diag = 0;
    double diag = 0;

    for(unsigned int i=0;i<M.n_rows;i++)
    {
        for(unsigned int j=0;j<M.n_cols;j++)
        {
            double v=0;
            for(unsigned int k=2; k<=M(i,j); k++)
                v+=log(k);
            with_diag+=v;
            if(i==j)
                diag+=v;
        }
    }

    without_diag = with_diag-diag;
}

inline
vec vech(const mat & M)
{
    unsigned int n=M.n_rows;
    vec v(n*(n+1)/2);
    
    unsigned int k=0;
    for(unsigned int i=0; i<n; i++)
    {
        for(unsigned int j=i;j<n; j++)
        {
            v(k++)=M(i,j);
        }
    }
    return(v);
}

inline 
mat unvech(const vec & v)
{
    unsigned int n=((sqrt(1+8*v.n_elem)+1)-1)/2;
    mat M(n,n);
    
    unsigned int k=0;
    for(unsigned int i=0; i<n; i++)
    {
        for(unsigned int j=i;j<n; j++)
        {
            M(i,j)=M(j,i)=v(k++);
        }
    }

    return(M);
}

inline
mat logit(mat & x)
{
    return(log(x)-log(1-x));
}

inline
mat sigmo(mat & x)
{
    return(1/(1+exp(-x)));
}

inline
double logit(double & x)
{
    return(log(x)-log(1-x));
}

inline
double sigmo(double & x)
{
    return(1/(1+exp(-x)));
}

inline
mat compute_B(colvec & beta, cube & covariates)
{
    // B definied as \sum \beta_k C_{:,:,k}
    mat B=zeros<mat>(covariates.n_rows,covariates.n_cols);

    for(unsigned int k=0;k<covariates.n_slices;k++)
    {
        B+=beta(k)*covariates.slice(k);
    }
    return B;
}

