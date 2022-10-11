Pkg.activate(".")
using LTEI_TA
using Test
using BenchmarkTools
using Dates
using Combinatorics
using MATLAB
mat"addpath('/Users/sbadredd/Desktop/low-rank-two-electron-integrals/SOLVER/chebfun-chebfun-c07b658')"

##
#Test time computation of Chebyshev coefficients
@testset "Test_coeff" begin
    file_name1="data/Molecules/CO2.xlsx"
    @time Atoms,Mol=mol_prop(file_name1)
    Nb=Mol.Nb;

    boite=[-10,10]
    prec=1e4 #Precision/Accuracy Nq2
    mu=0.1#Mu parameter
    Accuracy=3 #Nodes and weights of 1st Gaussian qudrature rule Nq1:
    t, w=lgwt(Accuracy,0,mu) #Nodes and weights of 1st Gaussian qudrature rule :
    T=lgwt(Int(prec),boite[1],boite[2]);
    @info "$(now()) Extract number of interpolation points"
    n1=number_INT(boite,t,100,Timelimit=1.,flag=false)

    ##Normal coeff
    @info "$(now()) Evaluate Coeff"
    @btime Coeff=create_coeff($t,$boite,$n1)


    ##Adaptatif coeff
    domain=hcat([-100:100]...);
    distance=36;
    tau=1e-10;
    prec=1e4 #Precision/Accuracy Nq2
    mu=0.1 #Mu parameter
    Accuracy=3
    r,r1=Pairs(Nb,Atoms)
    cluster=adaptatif(Atoms,r1,domain,distance,tau);
    #Nodes and weights of 1st Gaussian qudrature rule Nq1:
    t, w=lgwt(Accuracy,0,mu) #Nodes and weights of 1st Gaussian qudrature rule :
    n1=number_INT(cluster,t,100,Timelimit=5.);
    b_inf=zeros(1,length(cluster))
    b_sup=zeros(1,length(cluster))
    @inbounds for j in 1:length(cluster)
        b_inf[j]=cluster[j].boite[1]
        b_sup[j]=cluster[j].boite[2]
    end
    T=lgwt(Int(prec),b_inf,b_sup);
    m=collect(with_replacement_combinations(1:length(cluster),2));
    m= collect(hcat(m...)');
    @info "$(now()) Evaluate Coeff"
    @btime Coeff=create_coeff($t,$cluster,$n1,$m)
end
