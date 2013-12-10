
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
