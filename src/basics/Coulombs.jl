#License is MIT: https://github.com/sbadred/LTEI_TA.jl/blob/99b988ec2d84266e51a5a9b6a5acaf190c26e019/LICENSE

function lookfor(a::Array{Int64,2},b::Array{Int64,2})
        output= zeros(Int64,0)
        for j in 1:size(a,1)
                l=findall(vec(all(in.(a[j,:]',b),dims=2)) .==1);
                append!(output, l)
        end
        output
end

function Coulombs!(M,Mo,Nb,r1,r;tol=1e-12)
    rsim=[1:Nb 1:Nb]
    RowIdx1=lookfor(rsim,r1);
    r2=r1[setdiff(1:end, RowIdx1), :];
    RowIdx3=lookfor(r1,collect(r))
    RowIdx4=lookfor(r2,r1)
    rowrow=lookfor(r2,collect(r))
    y1=zeros(size(M,1),size(Mo,2));y2=similar(y1);
    mul!(y1,M,Mo[RowIdx3,:]);mul!(y2,M[:,RowIdx4],Mo[rowrow,:])
    y1=y1+y2
    y1[abs.(collect(y1)) .< tol] .=0
    return y1
end
