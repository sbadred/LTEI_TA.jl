#License is MIT: https://github.com/sbadred/LTEI_TA.jl/blob/99b988ec2d84266e51a5a9b6a5acaf190c26e019/LICENSE

"""This functions returns a set of different clusters of pairs of orbitals
according to a threshold """
##
mutable struct clus
    boite::Array{Int64,2}
    ao
    Percentage
end

function adaptatif(Atoms::Array{atom,1},r1::Array{Int64,2},domain::Array{Int64,2},distance::Number,tau::Float64)
    cluster=Array{clus, 1}()
    si=0;
 for i in 1:size(r1,1)
        l1,o1=search_(r1[i,1],Atoms);
        l2,o2=search_(r1[i,2],Atoms);
        alpha=Atoms[l1].Orbits[o1].expo[length(Atoms[l1].Orbits[o1].expo)]
        beta= Atoms[l2].Orbits[o2].expo[length(Atoms[l2].Orbits[o2].expo)]
        Rp=(alpha.*Atoms[l1].Geo+beta.*Atoms[l2].Geo)./(alpha+beta);
        dis_square=sum((Atoms[l1].Geo .- Atoms[l2].Geo).^2);
        cst=(alpha*beta)/(alpha+beta);
        if exp(-cst*(dis_square))>tau
            si+=1;
            Hirar1=(alpha+beta).*(Rp[1].-domain).^2;
            Hirar2=(alpha+beta).*(Rp[2].-domain).^2;
            Hirar3=(alpha+beta).*(Rp[3].-domain).^2;
            index1=findall(x -> x<distance, abs.(Hirar1))
            index2=findall(x -> x<distance, abs.(Hirar2))
            index3=findall(x -> x<distance, abs.(Hirar3))
            j_max=maximum([domain[index1[end]] domain[index2[end]] domain[index3[end]]])
            j_min=minimum([domain[index1[1]] domain[index2[1]] domain[index3[1]]])
            #boite=[-10 10]
            if abs(j_min)>j_max
                boite=[j_min abs(j_min)]
            else
                boite=[-abs(j_max) j_max];
            end
            if ~isempty(boite)
                #isVectorclus=findall(x -> x==1,vec(all(in.(boite,v),dims=2)))
                isVectorclus=[i for i=1:length(cluster) if isequal(cluster[i].boite,boite) ]
                if ~isempty(isVectorclus)
                    cluster[isVectorclus][1].ao=[cluster[isVectorclus][1].ao;r1[i,:]'];
                    cluster[isVectorclus][1].Percentage=cluster[isVectorclus][1].Percentage+(1*100)/length(r1);
                else
                    push!(cluster,clus(boite,r1[i,:]',(1*100)/length(r1)));
                end
            end
        end
    end
    return cluster
end

function optimal_box(cluster::Array{clus,1})
    b_inf=zeros(1,length(cluster))
    b_sup=zeros(1,length(cluster))
    @inbounds for j in 1:length(cluster)
           b_inf[j]=cluster[j].boite[1]
           b_sup[j]=cluster[j].boite[2]
    end
    boite=[-Int(maximum(b_sup)),Int(maximum(b_sup))]
    return boite
end
