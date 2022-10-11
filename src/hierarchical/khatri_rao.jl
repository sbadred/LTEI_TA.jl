#License is MIT: https://github.com/sbadred/LTEI_TA.jl/blob/99b988ec2d84266e51a5a9b6a5acaf190c26e019/LICENSE

"""function to compute khatri-rao product"""

function dotkron(B::Matrix{Float64},C::Matrix{Float64},D::Matrix{Float64})
    N = size(B,1)
    A=zeros(N,size(B,2)*size(C,2)*size(D,2));
    @inbounds Threads.@threads for n = 1:N
        @views A[n,:]=kron(B[n,:],C[n,:],D[n,:])
    end
    return A
end
function dotkron2(B::Matrix{Float64},C::Matrix{Float64},D::Matrix{Float64})
    N = size(B,2)
    A=zeros(size(B,1)*size(C,1)*size(D,1),N);
    @inbounds Threads.@threads for n = 1:N
        @views A[:,n]=kron(B[:,n],C[:,n],D[:,n])
    end
    return A
end
