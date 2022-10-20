#License is MIT: https://github.com/sbadred/LTEI_TA.jl/blob/99b988ec2d84266e51a5a9b6a5acaf190c26e019/LICENSE

#Evaluation elemetn-wise fourth-order tensor
using Statistics
include("../basics/search_.jl")
include("../basics/transform_.jl")
include("../basics/Toolbox.jl")
"""
##############################################################
solver_element_wise: Evaluation elemetn-wise fourth-order tensor
                     using function solver_element_wise(args...)
Inputs:
-Atoms: molecular data
-analytics: reference values
-vec_index: element to evaluate
-boite: computational box
-n1:  interpolation points

Output:
-Integral's values
##############################################################
"""
##
function g_construction(vec_index::Array{Int64,1},Atoms::Array{atom,1},T1::Array{Float64,2},t1::Array{Float64,2},t2::Array{Float64,1})
    l1,o1=search_(vec_index[1],Atoms);
    l2,o2=search_(vec_index[2],Atoms);
    coeff1,expo1,R1,const1=transform_(Atoms[l1].Orbits[o1].expo,Atoms[l1].Orbits[o1].coeff,Atoms[l1].Geo,Atoms[l2].Orbits[o2].expo,Atoms[l2].Orbits[o2].coeff,Atoms[l2].Geo);
    #Constant from the gaussian product rule
    constante1=broadcast(exp,-const1*( (Atoms[l1].Geo[1]-Atoms[l2].Geo[1])^2+(Atoms[l1].Geo[2]-Atoms[l2].Geo[2])^2+(Atoms[l1].Geo[3]-Atoms[l2].Geo[3])^2));

    #construcntion g_munu in each direction
    v1=zeros(size(T1,1),length(coeff1));v2=similar(v1);v3=similar(v1);
    A=similar(coeff1)
    broadcast!(*,A,coeff1,constante1)
    #First dimension x1
    AA=((t1 .- Atoms[l1].Geo[1]).^(Atoms[l1].Orbits[o1].couche[1]).*((t1 .-Atoms[l2].Geo[1]).^(Atoms[l2].Orbits[o2].couche[1])));
    B=zeros(length(A),length(t1))
    broadcast!(exp,B, -expo1.*(t1 .-R1[:,1]).^2);
    broadcast!(*,B,A,AA,B)
    mul!(v1,T1,broadcast!(*,B',t2,B'))
    #Second dimension x2
    AA=((t1 .-Atoms[l1].Geo[2]).^(Atoms[l1].Orbits[o1].couche[2])).*((t1 .-Atoms[l2].Geo[2]).^(Atoms[l2].Orbits[o2].couche[2]))
    broadcast!(exp,B, -expo1.*(t1 .-R1[:,2]).^2);
    broadcast!(*,B,AA,B)
    mul!(v2,T1,broadcast!(*,B',t2,B'))
    #third dimension x3
    AA=((t1 .-Atoms[l1].Geo[3]).^(Atoms[l1].Orbits[o1].couche[3])).*((t1 .-Atoms[l2].Geo[3]).^(Atoms[l2].Orbits[o2].couche[3]))
    broadcast!(exp,B, -expo1.*(t1 .-R1[:,3]).^2);
    broadcast!(*,B,AA,B)
    mul!(v3,T1,broadcast!(*,B',t2,B'))
    return v1,v2,v3
end

function solver_element_wise(Atoms::Array{atom,1},analytics::Array{Float64,1},vec_index::Array{Int64,2},n1::Array{Int64,1},T::Tuple{Array{Float64,1},Array{Float64,1}},boite::Array{Int64,1}
    ,Coeff::Array{Any,2},w::Array{Float64,1})
    Error=zeros(Float64,size(vec_index,1));
    maxi=Int(maximum(n1))
    t1=collect(T[1]');t2=T[2];
    T1=chebypoly(T[1],boite[1],boite[2],maxi);
    value=0
    Threads.@threads for index=1:size(vec_index,1)
        #Mu_nu
        v_munu_1,v_munu_2,v_munu_3=g_construction(vec_index[index,1:2],Atoms,T1,t1,t2)
        #kappa_lambda
        v_kl_1,v_kl_2,v_kl_3=g_construction(vec_index[index,3:4],Atoms,T1,t1,t2)
        L=zeros(Float64,length(w));
   #println("measure time")
   #@time begin
        Threads.@threads for i=1:length(w)
            F=zeros(Float64,size(v_munu_1,2),size(v_kl_1,2));
            A=collect((Coeff[i,2]*v_munu_1[1:n1[i],:])')*(collect(Coeff[i,1]')*v_kl_1[1:n1[i],:])
            B=collect((Coeff[i,2]*v_munu_2[1:n1[i],:])')*(collect(Coeff[i,1]')*v_kl_2[1:n1[i],:])
            C=collect((Coeff[i,2]*v_munu_3[1:n1[i],:])')*(collect(Coeff[i,1]')*v_kl_3[1:n1[i],:])
            broadcast!(*,F,A,B,C)
            L[i]=w[i]*sum(sum(F))
        end
    #end
        value=(2/sqrt(pi))*sum(L);
        Error[index]=abs(value-analytics[index])/abs(analytics[index])

    end
    return mean(Error),value


end
