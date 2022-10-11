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
