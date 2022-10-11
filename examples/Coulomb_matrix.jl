"""The aim of this program is to accelerate the tensor contractions
using  the long-range two-electron integrals tensor, the resutl of these
contractions will be the Coulomb matrix in the MO basis.
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
#"""In this example we will try to evaluate the coulomb matrix in molecular basis"""
using LTEI_TA
using Test
using Dates
using Combinatorics
using Printf
using MATLAB
using LinearAlgebra
mat"addpath('/Users/sbadredd/Desktop/low-rank-two-electron-integrals/SOLVER/chebfun-chebfun-c07b658')"

##
@testset "Coulomb Matrix " begin

    function main(file_name2::String,file_name3::String,mu::Float64,
        boite::Array{Int64,1},prec::Float64,
        Accuracy::Int64,Nb,r,r1)

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
        @timev C,si=create_M(Nb,T,Atoms,prec,n1,r1,r,Mo,boite,"false")

        ##Create the coulomb matrix in the molecular basis set
        #test
        @info "$(now()) Application of the compressed Tensor"
        #With  coeff compression
        @timev Coul=create_TEI(w,Coeff,C,n1)

        ##Relative error

        Matrix_Coul=extract_MO(file_name3,size(Mo,2))
        Error=norm(Coul-Matrix_Coul)/norm(Matrix_Coul)
        @printf("Error: %f  max: %d  ",Error,maximum(n1))        #end
        @test isapprox(Error,1e-4, atol=1e-4)
    end

    file_name1="data/Molecules/1-glycine.xlsx"
    @info "$(now()) Collecting  Data"
    @time Atoms,Mol=mol_prop(file_name1)

    Nb=Mol.Nb
    domain=hcat([-100:100]...);
    distance=15;
    tau=1e-10;
    r,r1=Pairs(Nb,Atoms)
    cluster=adaptatif(Atoms,r1,domain,distance,tau);
    boite=optimal_box(cluster)

    prec=1e4 #Precision/Accuracy Nq2
    mu=0.4#Mu parameter
    Accuracy=8 #Nodes and weights of 1st Gaussian qudrature rule Nq1:

    file_name2="data/Mo_coeff/1-glycine_mo.txt";
    file_name3="data/Coul_exchange/glycine1_0_4.ezfio.coulomb_MO"
    main(file_name2,file_name3,mu,boite,prec,Accuracy,Nb,r,r1)
end
