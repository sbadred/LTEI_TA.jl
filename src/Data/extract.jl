#License is MIT: https://github.com/sbadred/LTEI_TA.jl/blob/99b988ec2d84266e51a5a9b6a5acaf190c26e019/LICENSE

using Combinatorics
"""
##############################################################
read_TEI: creates the TEI of size (Nb^2 x Nb^2) tensor with
          analytical values from the fcidump file
Inputs:
-fcidump_filename
-orbital_count
##############################################################
"""

function read_TEI(fcidump_filename::String,orbital_count::Int)::Array{Float64,2}

    lines = readlines(fcidump_filename)


    two_electron_integral_tensor =
    zeros(orbital_count, orbital_count, orbital_count, orbital_count)

    Threads.@threads for line in lines[1:end]

        tokens = split(line)
        value = parse(Float64, tokens[1])
        i = parse(Int, tokens[2])
        j = parse(Int, tokens[3])
        k = parse(Int, tokens[4])
        l = parse(Int, tokens[5])


        for indices in Set([(i,j,k,l), (j,i,k,l), (i,j,l,k), (j,i,l,k),
            (k,l,i,j), (k,l,j,i), (l,k,i,j), (l,k,j,i)])

            two_electron_integral_tensor[indices...] = value
        end
    end



    return  collect(sparse(reshape(two_electron_integral_tensor,orbital_count^2,orbital_count^2)))
end

"""
##############################################################
extract: exract values from the fcidump file
Inputs:
-file_name
-Seed number
-Nb number of basis functions
##############################################################
"""
function extract(file_name::String,Seed::Int64,Nb::Int64)
    r=rand(1:Nb,Seed)
    lines=readlines(file_name)[r]
    analytics=zeros(Float64,Seed);
    vec_index=zeros(Int64,Seed,4);
    for p=1:Seed
        tokens = split(lines[p])
        value = parse(Float64, tokens[1])
        i = parse(Int, tokens[2])
        j = parse(Int, tokens[3])
        k = parse(Int, tokens[4])
        l = parse(Int, tokens[5])
        analytics[p]=value;
        vec_index[p,:]=[i j k l]
    end
    return analytics, vec_index
end

"""
##############################################################
extract_MO: Gives the matrix that contains the molecular
            orbitals coefficients size (Norb x Norb)
Inputs:
-file_name
-Norb : Number of molecular orbitals
##############################################################
"""
function extract_MO(file_name::String,Norb::Int64)
    Mo=zeros(Norb,Norb)
    lines=readlines(file_name)
    Threads.@threads for p=1:length(lines)
        tokens = split(lines[p])
        value = parse(Float64, tokens[1])
        i = parse(Int, tokens[2])
        j = parse(Int, tokens[3])
        Mo[i,j]=value
    end
    return(Mo)

end


"""
##############################################################
Pairs: gives the set of pairs of basis functions (μ,ν)
Input:
-Nb
-Atoms
-tau
##############################################################
"""
#optimal number of pairs of basis
function Pairs(Nb::Int64,Atoms::Array{LTEI_TA.atom,1};tau=1e-10)
    ##Create pairs of indices
    all_perm(xs, n) = vec(map(collect, Iterators.product(ntuple(_ -> xs, n)...)));
    r=collect(all_perm(1:Nb,2));
    r=collect(hcat(r...)');

    #Optimised r where we  take into consideration the symmetries
    r1=collect(with_replacement_combinations(1:Nb,2));
    r1=collect(hcat(r1...)');
    lis_r=[];
    for i in 1:size(r1,1)
       l1,o1=search_(r1[i,1],Atoms);
       l2,o2=search_(r1[i,2],Atoms);
       alpha=Atoms[l1].Orbits[o1].expo[length(Atoms[l1].Orbits[o1].expo)]
       beta= Atoms[l2].Orbits[o2].expo[length(Atoms[l2].Orbits[o2].expo)]
       Rp=(alpha.*Atoms[l1].Geo+beta.*Atoms[l2].Geo)./(alpha+beta);
       dis_square=sum((Atoms[l1].Geo .- Atoms[l2].Geo).^2);
       cst=(alpha*beta)/(alpha+beta);
       if exp(-cst*(dis_square))>tau
          push!(lis_r,r1[i,:])
       end
    end
    r1=collect(hcat(lis_r...)');
    return r,r1;
end
