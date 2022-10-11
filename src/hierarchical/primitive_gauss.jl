"""function to return the maximum number of linear combination of
primitive gaussian"""
##
function extract_prim(cluster::Array{clus,1},Atoms::Array{atom,1})
m=zeros(length(cluster))
for i=1:length(cluster)
        for k=1:size(cluster[i].ao,1)
                l1,o1=search_(cluster[i].ao[k,1],Atoms);
                l2,o2=search_(cluster[i].ao[k,2],Atoms);
                L=length(Atoms[l1].Orbits[o1].coeff)*length(Atoms[l2].Orbits[o2].coeff);
                if m[i]<L
                        m[i]=L;
                end
        end
end
return m
end
function extract_prim(cluster::Array{Int64,2},Atoms::Array{atom,1})
        m=0
        for k=1:size(cluster,1)
                l1,o1=search_(cluster[k,1],Atoms);
                l2,o2=search_(cluster[k,2],Atoms);
                L=length(Atoms[l1].Orbits[o1].coeff)*length(Atoms[l2].Orbits[o2].coeff);
                if m<L
                        m=L;
                end
        end

        return m
end
