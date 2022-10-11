"""
***************************************************************************************
Matt J (2022). Efficient Object-Oriented Kronecker Product
Manipulation (https://www.mathworks.com/matlabcentral/fileexchange/25969-efficient-object-oriented-kronecker-product-manipulation),
MATLAB Central File Exchange. Retrieved February 17, 2022
(some functions were converted to Julia language)
***************************************************************************************
"""


using Parameters
using Kronecker

function issquare(a)
#tests if argument is a square array (of whatever type).
if (size(a,1)==size(a,2)) & (length(size(a))<3)
 out=1;
else
 out=0;
end
out
end

function is1x1(a)
#tests if argument is a 1 x 1 object
if (size(a,1)==1) && issquare(a)==1
 out=1;
else
 out=0;
end
out=convert(Bool,out);
out
end
function isnum(a)
#Tests if argument is a numeric 1x1 object
#
#bool=isnum(a)
#
#Note: isnum(true)=true, isnum(NaN)=true
bool=is1x1(a) & (isa(a, Number) | isa(a, BitArray));
bool
end
function equaldims(size1,size2)
size1=size1'; size2=size2';
ii=findlast(!iszero, size1 .!==1);
jj=findlast(!iszero, size2 .!==1);
bool=!isempty(ii[1]==jj[1]) && (ii[1]==jj[1])&& all(size1[1:ii[1]]==size2[1:jj[1]]);
bool
end

function procShape(args...)

if length(args)<3;args[3]=0;end

   p=prod(args[1]);
   layers=prod(args[2])/p;

   if !convert(Bool,args[3])
    #z=findall(args[2], args[1])
    z=args[2]==args[1];

    if !z && z[1]==1;
      columnized=false;
   elseif args[2][1]==p;
      columnized=true;
    else
       print("Cannot determine reshaping rule")
    end

else #Transposed reshaping

    if args[2][1]==layers && args[2][2]==p && length(args[2])==2
      columnized=true;

   elseif equaldims(collect(args[2]),args[1])#implies 1 layer

        columnized=false;

    elseif equaldims( args[2][2:end], args[1])

        columnized=false;

    else

        print("Cannot determine reshaping rule")

    end

end
return columnized,Int(layers)
end

import Base.transpose
function transpose(M)
    Threads.@threads for ii in 1:M.numops
      M.opset[ii]=copy(M.opset[ii]') ;
   end
tmp=M.rangesizes;
M.rangesizes=M.domainsizes;
M.domainsizes=tmp;
M
end
function  full(M)
   if length(M.opinds)==2
      out=M.opset[1]⊗ M.opset[2];
   elseif length(M.opinds)==3
      out=M.opset[1]⊗M.opset[2]⊗M.opset[3];
end
end
import Base.size
function size(M)
   P=prod(M.rangesizes);Q=prod(M.domainsizes);
   return P,Q
end

function times(M,X)#for only scalarmask
if !isa(M,Kroneck) && isnum(M)
   X.scalarcoeff=X.scalarcoeff.*M;
   return X;

elseif !isa(X,Kroneck) && isnum(X) #implies M is KronProd
  tmp=X;
  X=M;
  X.scalarcoeff=X.scalarcoeff.*tmp;
  return X;
end
end
##
@with_kw mutable struct Kroneck
          opset :: Array{Any} =[]
          opinds :: Array{Any}=[]
          numops::Int =0
          eyemask=[]
          domainset=[]
          maxdim::Int =0
          scalarcoeff::Number =1
          scalarcumprods =1
          scalarmask=[]
          domainsizes=[]
          rangesizes=[]

   end
function write(M::Kroneck,args...)
   if length(args) <1;return;end
   if !isa(args[1],Array);args[1]=[args[1]];end
   if length(args)<2;opinds=1:length(args[1]);domainsizes=[];end
   if length(args)==3;opinds=args[2];domainsizes=[];end
   if length(args)>=4;opinds=args[2];domainsizes=args[3];M.scalarcoeff=args[4];end
   ##Weed out unused operands - stray operands will impact upon efficiency
   qq=opinds;ii=opinds';jj=opinds';
   opset=[qq];
   opinds=jj;
   #End Weed
   M.opset=args[1];
   M.opinds=opinds';
   M.numops=length(M.opset);
   M.maxdim=length(M.opinds);
   M.scalarcumprods=Int.(ones(size(M.opinds)));
   M.scalarmask=BitArray(zeros(1,M.maxdim));
   #M.domainsizes=domainsizes(:).';
   M.rangesizes=Int.(ones(size(M.domainsizes)));
   nColMap=Array{Int, 1}();
   for i in 1:length(M.opset)
     nColMap=push!(nColMap,size(M.opset[i],2));
   end
   nonfinMap=map(!,isfinite.(domainsizes));
   domainsizes=domainsizes';
   if (all(nonfinMap) && length(domainsizes)==1) || isempty(domainsizes)
     domainsizes=nColMap;
   else
     domainsizes(nonfinMap)=nColMap(nonfinMap);
   end
   M.domainsizes=domainsizes;

   if length(M.domainsizes)!=M.maxdim
     print("Inconsistency between domain sizes and operator index set.");
   end

   cc=0;
   for ii in M.opinds
     cc=cc+1;
     sz=size(args[1][ii]);
     if isnum(args[1][ii])
        M.rangesizes[cc]=M.domainsizes[cc];
        M.scalarmask[cc]=BitArray(1);
        M.scalarcumprods[cc]=M.opset[ii];
     elseif sz[2]==M.domainsizes[cc]
        push!(M.rangesizes,sz[1])
        #M.rangesizes[cc]=sz[1];
     else
        print(["DOMAINSIZE(' num2str(cc) ') not consistent with operand dimensions."]);
     end
   end
   M.scalarcumprods=cumprod(M.scalarcumprods);
   M.eyemask=similar(BitArray, (1,M.numops));
   M.domainset= nColMap;

   if any(x->x==1,M.scalarmask)
     U=[ M.opinds[M.scalarmask] ; M.domainsizes[M.scalarmask] ];
     U=unique(U',dims=1);
     opinds_shrink=unique(U[:,1]);

     if length(opinds_shrink)<length(U[:,1])
        print("A scalar member of OPSET is being used to represent two or more differently sized matrices")
     else
        M.eyemask[U[:,1]]=true;
        M.domainset[U[:,1]]=U[:,2]';
     end
   end
end
