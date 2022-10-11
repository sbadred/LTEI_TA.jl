#License is MIT: https://github.com/sbadred/LTEI_TA.jl/blob/99b988ec2d84266e51a5a9b6a5acaf190c26e019/LICENSE

function search_(N::Int64,Atoms::Array{atom,1})
    l=1;res=0;
    size_=0;minus=0;
    while l<=size(Atoms,1)
        size_=size_+size(Atoms[l].Orbits,1)
        if size_<N
            minus=minus+size(Atoms[l].Orbits,1)
            l=l+1
        else
            res=N-minus
            break
        end
    end
    return l,res
end
