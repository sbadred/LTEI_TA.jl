"""The aim of this program is to evaluate  the two-electron integrals
elements(six-dimensional integrals) for the erf interaction (long-range interaction)
where mu values are small for any given molecule for different basis-sets
We need to provide:
-File.xlsx of the molecule
-Get molecular properties
-Set a computational box [-b,b]
-Evaluate the number of interpolation points
Output: Value of the six-dimensional integral
"""
##
using LTEI_TA
using Test
using Dates
using Combinatorics
using Printf
using MATLAB
mat"addpath('/Users/sbadredd/Desktop/low-rank-two-electron-integrals/SOLVER/chebfun-chebfun-c07b658')"

##
function main(file_name1::String,
    boite::Array{Int64,1},
    t::Array{Float64,1},
    prec::Float64,
    w::Array{Float64,1},
    vec_index::Array{Int64,2},
    analytics:: Array{Float64,1},
    Atoms::Array{LTEI_TA.atom,1},n1,T)

    @info "$(now()) Create Chebyshev coefficient matrices"
    @timev Coeff=create_coeff(t,boite,n1)
    @info "$(now())Return the wanted numerical values of six-dimension integrals"
    @timev Error,value=solver_element_wise(Atoms ,analytics,vec_index,n1,T,boite,Coeff,w);
    @printf("Error: %d  max: %d  ",Error,maximum(n1))
    @test Error ≈1e-4 atol=1e-4
end

@testset "Novel method " begin

    @info "$(now()) Collecting data"
    #stucture  contatining molecule's properties: choose molecule and basis set
    file_name1="data/Molecules/CO2.xlsx"
    Atoms,Mol=mol_prop(file_name1)
    #Number of basis function
    Nb=Mol.Nb
    r,r1=Pairs(Nb,Atoms)
    domain=hcat([-100:100]...);
    distance=20;
    tau=1e-10;
    cluster=adaptatif(Atoms,r1,domain,distance,tau);
    boite=optimal_box(cluster)
    ## Numerical evaluation
    # Default Parameters for long-ranged interaction mu=0.5
    prec=1e4 #Precision/Accuracy Nq2
    mu=0.5 #Mu parameter
    Accuracy=7#Nodes and weights of 1st Gaussian qudrature rule Nq1:
    t, w=lgwt(Accuracy,0,mu) #Nodes and weights of 1st Gaussian qudrature rule :
    T=lgwt(Int(prec),boite[1],boite[2]);
    n1=number_INT(boite,t,100,Timelimit=1.,flag=false)    #Extract analytical values and the corresponded indices (random seed)
    @info "$(now()) Collecting analytical values"
    Seed=1000
    file_name2="data/Analytics/co2.ezfio.FCIDUMP_ao_erf"
    analytics,vec_index=extract(file_name2,Seed,Nb)
    #Or specefic index
    #vec_index=[24          10          17          28]
    #analytics=[-0.111975246633493]

    main(file_name1,boite,t,prec,w,vec_index,analytics,Atoms,n1,T)
end
