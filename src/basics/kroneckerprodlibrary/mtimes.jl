include("KroneckerProd.jl")
using SparseArrays
using  LinearAlgebra

function mtimes(L,R)
    if !isa(L,Kroneck) && isnum(L)  #means K=c*M
      c=L;M=R;
      out=M;
      out.scalarcoeff=out.scalarcoeff*c;
      return out;
    elseif !isa(R,Kroneck) && isnum(R) #means K=M*c
      M=L; c=R;
      out=M;
      out.scalarcoeff=out.scalarcoeff*c
      return out;

    elseif isa(L,Kroneck)  &  isa(R,Kroneck) #Means M=M1*M2;
         M1=L; M2=R;
         if !isequal(M1.domainsizes,M2.rangesizes)
           print("Incompatible sizes");
         end
         domainsizes=M2.domainsizes;
         scalarcoeff=M1.scalarcoeff*M2.scalarcoeff;
         if isequal(M1.opinds, M2.opinds)
           opinds=M2.opinds; zzz=unique(opinds);
           opset=Vector{Any}(undef, length(zzz))
           for ii in zzz
               opset[ii]=M1.opset[ii]*M2.opset[ii];
           end

         else
           opinds=1:length(M2.domainsizes);
           opset=Vector{Any}(undef, size(M.opinds));
           for ii in opinds
              opset[ii]=M1.opset[M1.opinds[ii]]*M2.opset[M2.opinds[ii]];
           end
        end
        out=Kroneck(opset,opinds,domainsizes,scalarcoeff);
        return out;

    elseif !isa(L,Kroneck) && (all(isa.(L, Number)) |isa(L, BitArray)) # means Y=X*M

        X=L; M=R;

        wasrowshaped,layers=procShape(collect(M.rangesizes'),size(X),1);

        Xt=reshape(X,(layers,:))';
        out=(full(transpose(M))*Xt)';

        if !wasrowshaped
         sz=M.domainsizes;
         if layers>1
           sz=[layers,sz]
         end
         out=reshape(out,sz);
        end
        return out;

      elseif !(isa(L,Kroneck)  &  (all(isa.(R, Number)) | isa(R, BitArray)))
        print("Undefined mult. operation involving KronProd");
      end
##
M=L; X=R;
restoresingle=false;

if all(isa.(X, Float32))
 for ii=M.opinds

    if issparse(M.opset{ii})
     X=convert(Array{Float64,2},X);
     restoresingle=true;
     break;
    end

 end
end
##
numel_domain=prod(M.domainsizes);
numel_range=prod(M.rangesizes);
##
wascolumnized,layers=procShape(M.domainsizes,size(X),0);
##

nn=M.maxdim;
current_dims=M.domainsizes;
if isnum(current_dims)
  current_dims[2]=1;
end
scalar=M.scalarcumprods[nn]*M.scalarcoeff;
if scalar==0 #Computationally easy case

    X=X[:];
    X[numel_range]=0;
    X=X[1:numel_range];
    X[1:numel_range]=0;

else
  if numel_domain<=numel_range && scalar!=1
    X=scalar*X;
  end
##LOOP OVER DIMENSION
@inbounds  @fastmath for ii in 1:nn
  X=reshape(copy(X),(current_dims[ii],:));
  if M.scalarmask[ii] #current operator represented by scalar
     X=X'; #nothing else necessary
  elseif issparse(M.opset[M.opinds[ii]]) && issparse(X)
     X=(X' * M.opset[M.opinds[ii]]');  #Note: no tranpose necessary
  else
     y=zeros(size(M.opset[M.opinds[ii]],1),size(X,2));
     mul!(y,copy(M.opset[M.opinds[ii]]),X);
     X=y';
  end
end
##

if numel_domain>numel_range && scalar!=1
  X=scalar*X;
end
end
if layers>1
     X=reshape(X,(layers,:))';
     #X=convert(Array{Float64,2},X);
end
#@show typeof(X)
if wascolumnized
  X=reshape(copy(X),(numel_range,:));
else
  X=reshape(copy(X),Tuple(vcat(M.rangesizes, layers)))
end
if restoresingle
  X=convert(Array{Float32,2},X);
end
#@show typeof(X)
out=X;
return out
end
