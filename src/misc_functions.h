
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
