include("KroneckerProd.jl")
include("mtimes.jl")


@time begin
A=rand(43,43)
S=A⊗ A⊗ A;
X=rand(43^3,33)
Y=S*X;
end

@time begin
#write(S1,collect(L))
Y1=mtimes(S1,X)
end

@time Y=S*X;
@time Y1=mtimes(S1,X)
