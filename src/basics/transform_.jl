#License is MIT: https://github.com/sbadred/LTEI_TA.jl/blob/99b988ec2d84266e51a5a9b6a5acaf190c26e019/LICENSE

function transform_(expo::Array{Float64,1},coeff::Array{Float64,1},R::Array{Float64,2},expo_::Array{Float64,1},coeff_::Array{Float64,1},R_::Array{Float64,2})
    c=kron(coeff,coeff_);
    e=expo .+ expo_';
    e=vec(e');
    ee=kron(expo,expo_);
    constante=ee./e;

    M1=(expo*R);
    M2=(expo_*R_);
    Rp=zeros(size(coeff,1)*size(coeff_,1),size(R,2));
    for j=1:size(Rp,2)
        A=M1[:,j] .+M2[:,j]';
        Rp[:,j]=vec(A')./e;
    end
    return c,e,Rp,constante
end
