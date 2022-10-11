"""The aim of this program is to accelerate the tensor contractions
using  the long-range two-electron integrals tensor, the resutl of these
contraction will be the Exchange matrix in the MO basis.
-mu values are small for any given molecule for different basis-sets
We need to provide:
-File.xlsx of the molecule
-Mo: Molecular coefficient matrix
-Get molecular properties
-Set a computational box [-b,b]
-Evaluate the number of interpolation points
Output: Value of the tensor contraction and Time executions
"""

##
#"""In this example we will try to evaluate the Exchange matrix in atomic basis"
using LTEI_TA
using Test
using Dates
using Combinatorics
using Printf
using MATLAB
using LinearAlgebra
mat"addpath('/Users/sbadredd/Desktop/low-rank-two-electron-integrals/SOLVER/chebfun-chebfun-c07b658')"
##
@testset "Exchange matrices" begin
    function main(file_name2::String,
        file_name3::String,mu::Float64,
        boite::Array{Int64,1},prec::Float64,
        Accuracy::Int64,Nb::Int64,r,r1)

        ##Molecular coefficient matrix
        D,Mo=orbitalss(Nb,file_name2);

        ## Numerical evaluation
        @info "$(now()) Pre-computation"
        t, w=lgwt(Accuracy,0,mu) #Nodes and weights of 1st Gaussian qudrature rule :
        T=lgwt(Int(prec),boite[1],boite[2]);
        @info "$(now()) Extract number of interpolation points"
        #best truncation box

        n1=number_INT(boite,t,100,Timelimit=1.,flag=false)

        @info "$(now()) Evaluate Coeff"
        @timev Coeff=create_coeff(t,boite,n1,flag=true)


        @info "$(now()) Creation of M matrix"
        #false without compression, true with compression
        @timev M=create_M(Nb,T,Atoms,prec,n1,r,boite)

        ##Create the coulomb matrix in the molecular basis set
        #test
        @info "$(now()) Application of the compressed Tensor"

        max_N1=maximum(n1)
        #AO
        Ml=reshape(M,max_N1,max_N1,max_N1,Nb,Nb);
        @timev K=create_exchange_ao(w,Coeff,Ml,Nb,size(D,2),n1,max_N1,D)
        ##Relative error
        MAtrix_Ex=extract_MO(file_name3,Nb)
        Error=norm(K-MAtrix_Ex)./norm(MAtrix_Ex)
        @printf("Error: %f  max: %d  ",Error,maximum(n1))
        @test isapprox(Error,1e-4, atol=1e-4)
    end

    file_name1="data/Molecules/ch4_vdz.xlsx"
    Atoms,Mol=mol_prop(file_name1)
    Nb=Mol.Nb
    r,r1=Pairs(Mol.Nb,Atoms)

    prec=1e4 #Precision/Accuracy Nq2
    mu=0.1  #Mu parameter
    Accuracy=4 #Nodes and weights of 1st Gaussian qudrature rule Nq1:

    @info "$(now()) Create optimal computational box"
    domain=hcat([-100:100]...);
    distance=15;
    tau=1e-10;
    cluster=adaptatif(Atoms,r1,domain,distance,tau);
    boite=optimal_box(cluster);


    file_name2="data/Mo_coeff/ch4_mo.txt";
    file_name3="data/Coul_exchange/ch4_01.ezfio.exchange_ao_"

    main(file_name2,file_name3,mu,boite,prec,Accuracy,Nb,r,r1)
end
