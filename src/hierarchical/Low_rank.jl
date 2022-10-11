#License is MIT: https://github.com/sbadred/LTEI_TA.jl/blob/99b988ec2d84266e51a5a9b6a5acaf190c26e019/LICENSE

"""Perform svd on given matrix"""
using TSVD
##
function Low_rank(vector::Array{Float64,3})
    vector=reshape(vector,size(vector,1),size(vector,2)*size(vector,3));
    k=rank(vector)
    U, s, V = TSVD.tsvd(A, k);
    return V,diagm(s)*collect(U');
end
