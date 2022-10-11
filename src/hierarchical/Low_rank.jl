"""Perform svd on given matrix"""
using TSVD
##
function Low_rank(vector::Array{Float64,3})
    vector=reshape(vector,size(vector,1),size(vector,2)*size(vector,3));
    k=rank(vector)
    U, s, V = TSVD.tsvd(A, k);
    return V,diagm(s)*collect(U');
end
