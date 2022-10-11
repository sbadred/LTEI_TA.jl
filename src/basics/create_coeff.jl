"""Create Chebyshev's coefficients matrices
##############################################################
Inputs:
-boite: Computational box
-t : quadrature nodes
-n1: number of interpolation nodes

Output:
-Chebyshev's coefficients matrices
##############################################################
"""

using LinearAlgebra
include("../basics/Toolbox.jl")
##
function sv_trunc(s::Array{Float64},tol)
	if tol==0.0
		return s
	else
		d = length(s)
		i=0
		weight = 0.0
		norm2 =dot(s,s)
		while (i<d) && weight<tol*norm2
			weight+=s[d-i]^2
			i+=1
		end
		return s[1:(d-i+1)]
	end
end
#simple case method
function create_coeff(t::Array{Float64,1},boite::Array{Int64,1},
	n1::Array{Int64,1};flag::Bool=true)

	if flag
		coeff=Array{Any}(undef, length(t),2);
	else
		coeff=Array{Any}(undef, length(t),1);
	end
	for i in 1:length(t)
		#function to be approximated
		fun(x,y)=exp(-t[i]^2 .*((x .-y).^2));
		A=chebynodes_grid(boite[1],boite[2],Int(n1[i]),boite[1],boite[2],Int(n1[i]));
		fcheby1=fun.(A[1],A[2]);
		COFF=interpspec2D_FFT(Int(n1[i]),Int(n1[i]),fcheby1);
		if flag
			U,eig,V=svd(COFF);
			s_trunc=sv_trunc(eig,1e-13);
			#indices=findall(eig.>=1e-7);
			coeff[i,1]=U[:,1:length(s_trunc)];
			coeff[i,2]=diagm(eig[1:length(s_trunc)])*V[:,1:length(s_trunc)]';
		else
			coeff[i]=COFF
		end
	end
	return coeff
end



##Hierarchical distribution
function create_coeff(t::Array{Float64,1},cluster::Array{clus,1},n1::Array{Int64,2},m::Array{Int64,2})
    coeff=Array{Any}(undef, length(t),2,size(m,1));
    Threads.@threads for p in 1:size(m,1)
        Threads.@threads for i in 1:length(t)
            #function to be approximated
            fun(x,y)=exp(-t[i]^2 .*((x .-y).^2));
            A=chebynodes_grid(cluster[m[p,1]].boite[1],cluster[m[p,1]].boite[2],Int(n1[i,m[p,1]]),cluster[m[p,2]].boite[1],cluster[m[p,2]].boite[2],Int(n1[i,m[p,2]]));
            fcheby1=fun.(A[1],A[2]);
            COFF=interpspec2D_FFT(Int(n1[i,m[p,1]]),Int(n1[i,m[p,2]]),fcheby1);
            U,eig,V=svd(COFF);
            indices=findall(eig.<1e-8);
            if !isempty(indices)
              coeff[i,1,p]=U[:,1:indices[1]];coeff[i,2,p]=diagm(eig[1:indices[1]])*V[:,1:indices[1]]'
            else
              coeff[i,1,p]=U;
              coeff[i,2,p]=diagm(eig)*V';
          end
        end

        end
        coeff
end
