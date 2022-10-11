
using TensorOperations
 """This function is used to construct the two-electron integrals tensor"""
#import library
include("../basics/kroneckerprodlibrary/KroneckerProd.jl")
include("../basics/kroneckerprodlibrary/mtimes.jl")

##
function KronProdMat(coeff1::Array{Float64,2},coeff2::Array{Float64,2},
        ww::Float64,M1::Array{Float64,2})
        #Create kronck objects
        U=Kroneck();
        V=Kroneck();
        write(U,[coeff1,coeff1,coeff1]);
        times(ww,U);
        write(V,[coeff2,coeff2,coeff2])
        L11=mtimes(U,M1);
        L22=mtimes(transpose(V),M1)
        y=zeros(size(L11,2),size(L22,2));
        mul!(y,L11',L22)
        return y;
end

#Coeff without compression
function KronProdMat(coeff1::Array{Float64,2},ww::Float64,M1::Array{Float64,2})
        #Create kronck objects

        A=Kroneck();
        write(A,[coeff1,coeff1,coeff1]);
        times(ww,A);
        L11=mtimes(A,M1);
        y=zeros(size(L11,2),size(L11,2));

        mul!(y,M1',L11)
        return y;
end

function KronProdMat(coeff1::Array{Float64,2},coeff12::Array{Float64,2},coeff13::Array{Float64,2},coeff21::Array{Float64,2},coeff22::Array{Float64,2},coeff23::Array{Float64,2},ww::Float64,M1::Array{Float64,2},M2::Array{Float64,2})
        U=Kroneck();
        V=Kroneck();
        write(U,[coeff13,coeff12,coeff1]);
        times(ww,U);
        write(V,[coeff23,coeff22,coeff21]);
        L11=mtimes(transpose(U),collect(M1'));
        L22=mtimes(V,collect(M2'))
        y=zeros(size(L11,2),size(L22,2));
        mul!(y,L11',L22)
        return y;
end

function KronProdMat(coeff1::Array{Float64,2},coeff2::Array{Float64,2}
        ,ww::Float64,M1::Array{Float64,2},M2::Array{Float64,2})
#Create kronck objects
        U=Kroneck();
        V=Kroneck();
        write(U,[coeff1,coeff1,coeff1]);
        times(ww,U);
        write(V,[coeff2,coeff2,coeff2])
        L11=mtimes(U,M1);
        L22=mtimes(transpose(V),M2)
        y=zeros(size(L11,2),size(L22,2));
        mul!(y,L11',L22)
        return y;
end



##
#Without coeff compression
function create_TEI(w::Array{Float64,1},coeff::Array{Any,2},M::Array{Float64,4}
        ,n1::Array{Int64,1};flag::Bool=true)
        L=Array{Any}(undef, length(w));
        for i in 1:length(w)
                Mu=reshape(M[1:n1[i],1:n1[i],1:n1[i],1:size(M,4)],(n1[i]*n1[i]*n1[i],size(M,4)));
                L[i]=KronProdMat(coeff[i,1]',w[i],Mu);
        end

        B=2/sqrt(pi)*(sum(L));
        B
end

##

function create_TEI(w::Array{Float64,1},coeff::Array{Any,2},M::Array{Float64,4},n1::Array{Int64,1})
        L=Array{Any}(undef, length(w));
        Threads.@threads for i in 1:length(w)

                Mu=reshape(M[1:n1[i],1:n1[i],1:n1[i],1:size(M,4)],(n1[i]*n1[i]*n1[i],size(M,4)));
                L[i]=KronProdMat(copy(coeff[i,1]'),copy(coeff[i,2]'),w[i],convert(Array{Float64,2},Mu));
        end

        B=2/sqrt(pi)*(sum(L));
        B
end
##
#TEI_vec
function create_TEI(w::Array{Float64,1},coeff::Array{Any,2},
                    M2::Array{Float64,4},M1::Array{Float64,4},
                    n1::Array{Int64,1})
        L=Array{Any}(undef, length(w));
        Threads.@threads for i in 1:length(w)
                Mu2=reshape(M2[1:n1[i],1:n1[i],1:n1[i],1],(n1[i]*n1[i]*n1[i],1));
                Mu1=reshape(M1[1:n1[i],1:n1[i],1:n1[i],1:size(M1,4)],(n1[i]*n1[i]*n1[i],size(M1,4)));
                L[i]=KronProdMat(copy(coeff[i,1]'),copy(coeff[i,2]'),w[i],Mu1,Mu2);
        end
        B=2/sqrt(pi)*(sum(L));
        B
end
##Adaptive TEI
function create_TEI(w::Array{Float64,1},coeff::Array{Any,3},n1::Array{Int64,2},M::Array{Any,1},L_modif::Array{Any,2},m::Array{Int64,2})
        L=Array{Any}(undef, length(w),size(m,1));
        Threads.@threads for p in 1:size(m,1)
                Threads.@threads for i in 1:length(w)
                        coeff1=L_modif[1,m[p,1]][:,1:n1[i,m[p,1]]];coeff2=L_modif[2,m[p,1]][:,1:n1[i,m[p,1]]];coeff3=L_modif[3,m[p,1]][:,1:n1[i,m[p,1]]];
                        #right-hand multiplication
                        coeff1_T=L_modif[1,m[p,2]][:,1:n1[i,m[p,2]]]';coeff2_T=L_modif[2,m[p,2]][:,1:n1[i,m[p,2]]]';coeff3_T=L_modif[3,m[p,2]][:,1:n1[i,m[p,2]]]';

                       #Evaluation
                       B1=coeff1*coeff[i,2,p]';B2=coeff2*coeff[i,2,p]';B3=coeff3*coeff[i,2,p]';
                       B1_T=coeff[i,1,p]'*collect(coeff1_T);B2_T=coeff[i,1,p]'*collect(coeff2_T);B3_T=coeff[i,1,p]'*collect(coeff3_T);

                       L[i,p]=KronProdMat(B1,B2,B3,B1_T,B2_T,B3_T,w[i],M[m[p,1]],M[m[p,2]]);
               end
       end
       Res=Array{Any}(undef, size(m,1));
       Threads.@threads for j=1:size(m,1)
               if m[j,1]!=m[j,2]
                       Res[j]=2/sqrt(pi)*(sum(L[:,j])+sum(L[:,j]'));
               else
                       Res[j]=2/sqrt(pi)*(sum(L[:,j]));
               end
       end
       B=sum(Res);
       B
end
##
function create_exchange(w::Array{Float64,1},coeff::Array{Any,2},Ml::Array{Float64,2},
        N::Int64,n1::Array{Int64,1},max_N1::Int64,Mo::Array{Float64,2})
        K=zeros(N,N);
        for i=1:N
                for j=1:N
                        L=Ml*Mo[:,(j-1)*N+i]
                        L=reshape(L,max_N1,max_N1,max_N1,:)
                        #V=copy(Ml[:,:,:,(j-1)*N+i]);
                        K[i,j]=create_TEI(w,coeff,L,n1)[1]
                end
        end

        return K
end

function create_exchange(w::Array{Float64,1},coeff::Array{Any,2},Ml::Array{Float64,4},
        N::Int64,n1::Array{Int64,1},max_N1::Int64,Mo::Array{Float64,2})
        K=zeros(N^2);
        @tensoropt L[x,y,z,e]:=Ml[x,y,z,d] * Mo[d,e]
        #K=create_TEI(w,coeff,L,n1)
        Threads.@threads for i=1:N^2
                V=copy(L[:,:,:,i]);
                K[i]=create_TEI(w,coeff,V[:,:,:,:],n1)[1]
        end

return K
end

function create_exchange_ao(w::Array{Float64,1},coeff::Array{Any,2},Ml::Array{Float64,5},
        Nb::Int64,Norb::Int64,n1::Array{Int64,1},max_N1::Int64,D::Array{Float64,2})
        K=zeros(Nb,Nb);
        @tensoropt L[x,y,z,μ,Norb]:=Ml[x,y,z,μ,λ] * D[λ,Norb]

        #K=create_TEI(w,coeff,L,n1)
        Threads.@threads for j=1:Norb
                K .=K +create_TEI(w,coeff,L[:,:,:,:,j],n1)
        end
return K
end
