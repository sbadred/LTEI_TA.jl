function orbitalss(Nb::Int64,file_name2::String)
    Mo_read=readdlm(file_name2);
    Mo=zeros((Nb,Int(size(Mo_read,1)/Nb)))

    i=1;j=1;
    while j<=size(Mo,2)
    Mo[:,j]=Mo_read[i:Nb+i-1];
    i+=Nb;j+=1;
    end
    L_=zeros(Nb^2,size(Mo,2))
    for j in 1:size(Mo,2)
        L_[:,j]=kron(Mo[:,j],Mo[:,j]);
    end
    return Mo,L_
end
function orbitalss(Nb::Int64,file_name2::String)
    Mo_read=readdlm(file_name2);
    Mo=zeros((Nb,Int(size(Mo_read,1)/Nb)))

    i=1;j=1;
    while j<=size(Mo,2)
    Mo[:,j]=Mo_read[i:Nb+i-1];
    i+=Nb;j+=1;
    end
    L_=zeros(Nb^2,size(Mo,2))
    for j in 1:size(Mo,2)
        L_[:,j]=kron(Mo[:,j],Mo[:,j]);
    end
    return Mo,L_
end
