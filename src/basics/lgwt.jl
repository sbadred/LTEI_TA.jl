using FastGaussQuadrature
"""mapping gauss quadrature nodes weights from [-1 1] to [a b]:
##############################################################
Input:
-N: number of quadrature points
-a,b: the integrals bounds

Output:
-x: Quadrature nodes
-w: Quadrature weigths
##############################################################
"""
function lgwt(N::Int64,a::Number,b::Number)
    x,w=gausslegendre(N);
    w=((b .- a)./2).*w;
    x=(a.*(1 .-x)+b.*(1 .+x))/2;
    return x,w
end
function lgwt(N::Int64,a::Array{Float64,2},b::Array{Float64,2})
    x,w=gausslegendre(N);
    w=((b .- a)./2).*w;
    x=(a.*(1 .-x)+b.*(1 .+x))/2;
    return x,w
end
