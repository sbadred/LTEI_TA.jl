"""This function returns the matrix M in Nb^2xn^3 used for the evaluation of
the TEI tensor"""
##librairies
using TSVD
using Kronecker
using LinearAlgebra
using LowRankApprox
include("../basics/search_.jl")
include("../basics/transform_.jl")
include("../basics/Toolbox.jl")
include("../basics/solver_element_wise.jl")
include("../basics/Coulombs.jl")
include("../hierarchical/primitive_gauss.jl")
include("../hierarchical/khatri_rao.jl")
include("../hierarchical/Low.jl")
##
function create_M(Nb::Int64,T::Tuple{Array{Float64,1},Array{Float64,1}},
        Atoms::Array{atom,1},prec::Float64,n1::Array{Int64,1},
        r::Array{Int64,2},boite::Array{Int64,1},flag::Bool=false;low::Bool=false)

        maxi=Int(maximum(n1))
        t1=collect(T[1]');t2=T[2];
        T1=chebypoly(T[1],boite[1],boite[2],maxi);
        M=zeros(maxi*maxi*maxi,size(r,1));
        Threads.@threads for k=1:size(r,1)
                l1,o1=search_(r[k,1],Atoms);
                l2,o2=search_(r[k,2],Atoms);
                coeff1,expo1,R1,const1=transform_(Atoms[l1].Orbits[o1].expo,Atoms[l1].Orbits[o1].coeff,Atoms[l1].Geo,Atoms[l2].Orbits[o2].expo,Atoms[l2].Orbits[o2].coeff,Atoms[l2].Geo);

                #Constant from the gaussian product rule
                constante1=broadcast(exp,-const1*( (Atoms[l1].Geo[1]-Atoms[l2].Geo[1])^2+(Atoms[l1].Geo[2]-Atoms[l2].Geo[2])^2+(Atoms[l1].Geo[3]-Atoms[l2].Geo[3])^2));

                v1,v2,v3=g_construction(r[k,:],Atoms,T1,t1,t2)

                L=zeros(maxi*maxi*maxi,size(v1,2));
                for j=1:size(v1,2)
                        L[:,j]=kron(v1[:,j],kron(v2[:,j],v3[:,j]))
                end
                M[:,k]=sum(L, dims = 2);

        end
        if !low
                if !flag
                        M=reshape(M,maxi,maxi,maxi,size(r,1));
                        return(M)
                else
                        #interpolative decomposition
                        opts = LRAOptions(maxdet_tol=0., sketch_randn_niter=1)
                        opts.rtol = 1e-8
                        V=idfact(M, opts);
                        S=V[:sk];perm=V[:p];T_M=V[:T];
                        k=length(S)
                        U=reshape(M[:,S],maxi,maxi,maxi,length(S));
                        LL=[Matrix(I,k,k) T_M];
                        LL_=zeros(length(S),size(r,1));
                        Threads.@threads  for i=1:size(r,1)
                                idd=findall(perm .==i)
                                LL_[:,i]=LL[:,idd];
                        end



                        return U,LL_




                end
        else
                L_modif=Array{Any}(undef,3);
                m_modif=extract_prim(r,Atoms)
                maxi=Int(maximum(n1));
                T1=chebypoly(T[1],boite[1],boite[2],maxi)
                t1=collect(T[1]');t2=T[2];
                vec1=zeros(maxi,Int(m_modif),size(r,1));
                vec2=similar(vec1);
                vec3=similar(vec1);
                Threads.@threads for k=1:size(r,1)
                        v1,v2,v3=g_construction(r[k,:],Atoms,T1,t1,t2)
                        #Modif
                        vec1[:,:,k]=[v1 zeros(maxi,Int(m_modif)-size(v1,2))];
                        vec2[:,:,k]=[v2 zeros(maxi,Int(m_modif)-size(v2,2))];
                        vec3[:,:,k]=[v3 zeros(maxi,Int(m_modif)-size(v3,2))];
                end
                V1,L1,s1=Low_rank(vec1)
                L_modif[1]=L1;
                V2,L1,s2=Low_rank(vec2)
                L_modif[2]=L1;
                V3,L1,s3=Low_rank(vec3)
                L_modif[3]=L1;
                #Fill M
                M1_=zeros(size(V1,2)*size(V2,2)*size(V3,2),size(r,1));
                V1_=reshape(V1,Int(m_modif),size(r,1),size(V1,2));
                V2_=reshape(V2,Int(m_modif),size(r,1),size(V2,2));
                V3_=reshape(V3,Int(m_modif),size(r,1),size(V3,2));
                Threads.@threads for k=1:size(r,1)
                        L=zeros(size(V1,2)*size(V2,2)*size(V3,2),Int(m_modif));
                        for j=1:Int(m_modif)
                                L[:,j]=kron(V1_[j,k,:],kron(V2_[j,k,:],V3_[j,k,:]))

                        end
                        M1_[:,k]=sum(L,dims=2)
                end


                return M1_,L_modif
        end

end
#For coulomb
function create_M(Nb::Int64,T::Tuple{Array{Float64,1},Array{Float64,1}},
        Atoms::Array{atom,1},prec::Float64,n1::Array{Int64,1},
        r1::Array{Int64,2},r::Array{Int64,2},Mo,boite::Array{Int64,1},flag::String
        ;low::Bool=false)

        maxi=Int(maximum(n1))
        t1=collect(T[1]');t2=T[2];
        T1=chebypoly(T[1],boite[1],boite[2],maxi);
        M=zeros(maxi*maxi*maxi,size(r1,1));
        if !low
                Threads.@threads for k=1:size(r1,1)
                        l1,o1=search_(r1[k,1],Atoms);
                        l2,o2=search_(r1[k,2],Atoms);
                        coeff1,expo1,R1,const1=transform_(Atoms[l1].Orbits[o1].expo,Atoms[l1].Orbits[o1].coeff,Atoms[l1].Geo,Atoms[l2].Orbits[o2].expo,Atoms[l2].Orbits[o2].coeff,Atoms[l2].Geo);

                        #Constant from the gaussian product rule
                        constante1=broadcast(exp,-const1*( (Atoms[l1].Geo[1]-Atoms[l2].Geo[1])^2+(Atoms[l1].Geo[2]-Atoms[l2].Geo[2])^2+(Atoms[l1].Geo[3]-Atoms[l2].Geo[3])^2));

                        v1,v2,v3=g_construction(r1[k,:],Atoms,T1,t1,t2)

                        L=zeros(maxi*maxi*maxi,size(v1,2));
                        for j=1:size(v1,2)
                                L[:,j]=kron(v1[:,j],kron(v2[:,j],v3[:,j]))
                        end
                        M[:,k]=sum(L, dims = 2);
                end



                if flag =="false"
                        si=size(M);
                        M=Coulombs!(M,Mo,Nb,r1,r)
                        M=reshape(M,maxi,maxi,maxi,size(M,2))
                        return M,si;
                else
                        #interpolative decomposition
                        opts = LRAOptions(maxdet_tol=0., sketch_randn_niter=1)
                        opts.rtol = 1e-8
                        V=idfact(M, opts);
                        S=V[:sk];perm=V[:p];T_M=V[:T];
                        k=length(S)
                        U=reshape(M[:,S],maxi,maxi,maxi,length(S));
                        LL=[Matrix(I,k,k) T_M];
                        LL_=zeros(length(S),size(r,1));
                        Threads.@threads  for i=1:size(r,1)
                                idd=findall(perm .==i)
                                LL_[:,i]=LL[:,idd];
                        end
                end
                return U,LL_
        else  #khatri_rao
                L_modif=Array{Any}(undef,3);
                m_modif=extract_prim(r1,Atoms)
                maxi=Int(maximum(n1));
                T1=chebypoly(T[1],boite[1],boite[2],maxi)
                t1=collect(T[1]');t2=T[2];
                vec1=zeros(maxi,Int(m_modif),size(r1,1));
                vec2=similar(vec1);
                vec3=similar(vec1);
                Threads.@threads for k=1:size(r1,1)
                        v1,v2,v3=g_construction(r1[k,:],Atoms,T1,t1,t2)
                        #Modif
                        vec1[:,:,k]=[v1 zeros(maxi,Int(m_modif)-size(v1,2))];
                        vec2[:,:,k]=[v2 zeros(maxi,Int(m_modif)-size(v2,2))];
                        vec3[:,:,k]=[v3 zeros(maxi,Int(m_modif)-size(v3,2))];
                end
                V1,L1,s1=Low_rank(vec1)
                L_modif[1]=L1;
                V2,L1,s2=Low_rank(vec2)
                L_modif[2]=L1;
                V3,L1,s3=Low_rank(vec3)
                L_modif[3]=L1;
                #Fill M
                M1_=zeros(size(V1,2)*size(V2,2)*size(V3,2),size(r1,1));
                V1_=reshape(V1,Int(m_modif),size(r1,1),size(V1,2));
                V2_=reshape(V2,Int(m_modif),size(r1,1),size(V2,2));
                V3_=reshape(V3,Int(m_modif),size(r1,1),size(V3,2));
                Threads.@threads for k=1:size(r1,1)
                        L=zeros(size(V1,2)*size(V2,2)*size(V3,2),Int(m_modif));
                        for j=1:Int(m_modif)
                                L[:,j]=kron(V1_[j,k,:],kron(V2_[j,k,:],V3_[j,k,:]))

                        end
                        M1_[:,k]=sum(L,dims=2)
                end


                return M1_,L_modif
        end

end

##For coulomb with adaptaitf
function create_M(Nb::Int64,cluster::Array{clus,1},
        T::Tuple{Array{Float64,2},Array{Float64,2}},
        Atoms::Array{atom,1},precision::Float64,n1::Array{Int64,2},
        Mo::Array{Float64,2},r::Array{Int64,2};flag::Bool=true)
        #Modif
        @time begin
        M = Vector{Any}(undef, length(cluster))
        L_modif=Array{Any}(undef,3,length(cluster));
        m_modif=extract_prim(cluster,Atoms)
        Threads.@threads for i=1:length(cluster)
                maxi=Int(maximum(n1[:,i]));
                T1=chebypoly(T[1][:,i],cluster[i].boite[1],cluster[i].boite[2],maxi)
                t1=collect(T[1][:,i]');t2=T[2][:,i];
                vec1=zeros(maxi,Int(m_modif[i]),size(cluster[i].ao,1));
                vec2=similar(vec1);
                vec3=similar(vec1);
                Threads.@threads for k=1:size(cluster[i].ao,1)

                        v1,v2,v3=g_construction(cluster[i].ao[k,:],Atoms,T1,t1,t2)

                        #Modif
                        vec1[:,:,k]=[v1 zeros(maxi,Int(m_modif[i])-size(v1,2))];
                        vec2[:,:,k]=[v2 zeros(maxi,Int(m_modif[i])-size(v2,2))];
                        vec3[:,:,k]=[v3 zeros(maxi,Int(m_modif[i])-size(v3,2))];
                end

                #Modif


                #if (maxi<Int(m[i])*size(cluster[i].ao,1))
                        V1,L1=Low_rank(vec1)
                        L_modif[1,i]=L1;

                        V2,L1=Low_rank(vec2)
                        L_modif[2,i]=L1;

                        V3,L1=Low_rank(vec3)
                        L_modif[3,i]=L1;

              #=else
                        F=svd(reshape(vector1,maxi,Int(m_modif[i])*size(cluster[i].ao,1)),full=true);
                        indices=findall(F.S .<1e-8);
                        V1=collect(F.Vt[1:indices[end],:]');
                        L_modif[1,i]=F.U[:,1:indices[end]]*diagm(F.S[1:indices[end]]);
                        F=svd(reshape(vector2,maxi,Int(m_modif[i])*size(cluster[i].ao,1)));
                        indices=findall(F.S .<1e-8);
                        V2=collect(F.Vt[1:indices[end],:]');
                        L_modif[2,i]=F.U[:,1:indices[end]]*diagm(F.S[1:indices[end]]);
                        F=svd(reshape(vector3,maxi,Int(m_modif[i])*size(cluster[i].ao,1)));
                        indices=findall(F.S .<1e-8);
                        V3=collect(F.Vt[1:indices[end],:]');
                        L_modif[3,i]=F.U[:,1:indices[end]]*diagm(F.S[1:indices[end]]);
                end=#
                #Fill M
                M1_=zeros(size(V1,2)*size(V2,2)*size(V3,2),size(cluster[i].ao,1));
                V1_=reshape(V1,Int(m_modif[i]),size(cluster[i].ao,1),size(V1,2));
                V2_=reshape(V2,Int(m_modif[i]),size(cluster[i].ao,1),size(V2,2));
                V3_=reshape(V3,Int(m_modif[i]),size(cluster[i].ao,1),size(V3,2));

                Threads.@threads for k=1:size(cluster[i].ao,1)
                        L=zeros(size(V1,2)*size(V2,2)*size(V3,2),Int(m_modif[i]));
                        for j=1:Int(m_modif[i])
                                #v1=T1*((coeff1[j]*constante1[j]*((t1 .- Atoms[l1].Geo[1]).^(Atoms[l1].Orbits[o1].couche[1])).*((t1 .-Atoms[l2].Geo[1]).^(Atoms[l2].Orbits[o2].couche[1]))).*exp.(-expo1[j]*(t1 .-R1[j,1]).^2))'
                                #v2=T1*(((t1 .-Atoms[l1].Geo[2]).^(Atoms[l1].Orbits[o1].couche[2])).*((t1 .-Atoms[l2].Geo[2]).^(Atoms[l2].Orbits[o2].couche[2])).*exp.(-expo1[j]*(t1 .-R1[j,2]).^2))';
                                #v3=T1*(((t1 .-Atoms[l1].Geo[3]).^(Atoms[l1].Orbits[o1].
                                L[:,j]=kron(V1_[j,k,:],kron(V2_[j,k,:],V3_[j,k,:]))
                                #L[:,j]= collect(v1[:,j]'⊗v2[:,j]'⊗ v3[:,j]');
                        end
                        M1_[:,k]=sum(L,dims=2)

                end


                #=@time begin
                for p=1:Int(m_modif[i])
                        M1_+=dotkron(V1_[p,:,:],V2_[p,:,:,],V3_[p,:,:]);
                end
                end=#
                #si=size(M1_)
                if flag
                        y1=Coulombs!(M1_,collect(Mo),Nb,cluster[i].ao,collect(r))
                        M[i]=collect(y1');
                else
                        M[i]=M1_;
                end
                @show cluster[i].boite
        end
end
        return M,L_modif

end
