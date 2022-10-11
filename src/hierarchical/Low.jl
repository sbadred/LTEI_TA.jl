"""Perform svd on given matrix"""

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
function Low_rank(vector::Array{Float64,3})
    #=vector=reshape(vector,size(vector,1),size(vector,2)*size(vector,3));
    k=rank(vector)
    U, s, V = TSVD.tsvd(vector, k);=#
    F=svd(reshape(vector,size(vector,1),size(vector,2)*size(vector,3)))
    indices=findall(F.S .<=1e-8);
	if isempty(indices)
		return F.V,diagm(F.S)*collect(F.U'),F.S;
	else
		V1=F.V[:,1:indices[1]];
		L1=diagm(F.S[1:indices[1]])*collect(F.U[:,1:indices[1]]');
		return V1,L1,F.S[1:indices[1]];
	end
end
