using FFTW
"""Chebyshev nodes
# Given an interval (a,b), calculates the associated n Chebyshev nodes.
#
# Input:
# a,b= interval of interpolation : integers
# n=number of interpolation points: integer
#
# Output:
# nodes= nodes of chebyshev interpolation: vector in R^(n)"""
##--------------------------------------------------------------------------
function  chebynodes(a::Int64, b::Int64, n::Int64)

nodes = zeros(1,n);

for ii in 1:n
    #nodes[1,ii] = 0.5 * (b + a) - 0.5 * (a - b) * cos(((2*ii -1) * pi)/(2*(n)));
    nodes[1,ii] = 0.5 * (b + a) - 0.5 * (a - b) * cos(((ii - 1) * pi)/(n-1));

end
return nodes

end

## Chebyshev polynomials
"""Calculates the polynomials functions of N degree, interval [a,b]
# Input:
# x= vector
# a,b= integers interval of interpolation.
# N=number of interpolation points
#
# Output:
# chebyplys= chebyshev's polynomials evaluates on the vector x: matrix in
# R^(size(x,1)xN)"""
#--------------------------------------------------------------------------
function chebypoly(x::Array{Float64,1},a::Int64,b::Int64,N::Int64)
    # Necessary transformation to ensure that x lies in [-1,1]
    x = (2 .* (x .-a) / (b-a) ) .- 1;
    #Initialization
    chebypolys=zeros(N,length(x))
    #T_1(x)=1
    chebypolys[1,:] .=1
    #T_2(x)=x
    if N>1
        chebypolys[2,:] .=x;
    end
    # for n>2, T_n(x) = 2*x*T_{n-1}(x) - T_{n-2}(x)
    if N>2
        for k in 3:N
            chebypolys[k,:] .= 2 * x.* chebypolys[k-1,:] - chebypolys[k-2,:]
        end
    end

    chebypolys
end
function chebypoly(x::Number,a::Int64,b::Int64,N::Int64)
    # Necessary transformation to ensure that x lies in [-1,1]
    x = (2 .* (x .-a) / (b-a) ) .- 1;
    #Initialization
    chebypolys=zeros(N,length(x))
    #T_1(x)=1
    chebypolys[1,:] .=1
    #T_2(x)=x
    if N>1
        chebypolys[2,:] .=x;
    end
    # for n>2, T_n(x) = 2*x*T_{n-1}(x) - T_{n-2}(x)
    if N>2
        for k in 3:N
            chebypolys[k,:] .= 2 * x.* chebypolys[k-1,:] - chebypolys[k-2,:]
        end
    end

    chebypolys
end
## 2D Chebyshev grid
"""Calculates the [X,Y] grid obtained by the mesh of two vx, vy vectors.
#Input:
# a1,b1,a2,b2=intervals of 2D interpolation
# n1,n2=number of interpolation points for each dimension
#
#OUtput:
#[X,Y]= chebyshev's 2D grid: matrix in R^(n1,n2)"""
#--------------------------------------------------------------------------

function chebynodes_grid(a1::Int64, b1::Int64, n1::Int64, a2::Int64, b2::Int64, n2::Int64)
meshgrid(x,y) = (repeat(x',length(y),1),repeat(y,1,length(x)))
vx = chebynodes(a1,b1,n1);
vy = chebynodes(a2,b2,n2);

A= meshgrid(vx',vy');
return A
end

## 2D Spectral interpolation using FFT
"""n1,n2: number of Chebyshev nodes
# a1,a2: intervall [a1,b1], [a2,b2]
# b1,b2: intervall [a1,b1], [a2,b2]
# f: function to be interpolated, forme f
#
#Output:
# coefficients2= matrix containing chebyshev's coefficients in R^(n1xn2)
Reference:
"""

#--------------------------------------------------------------------------
"""Calculation of Chebyshev coefficients using the FFT method
***************************************************************************************
*    Title: SpectralInterpolation
*    Author: Edith Viau
*    Date: 2015
*    Availability: https://github.com/xuzhibo/SpectralInterpolation
*
***************************************************************************************
"""
function coeff_fft(value::Array{Float64,2})
    N=size(value,2)-1
    #coeff=real(fft([value';value[:,2:end-1]'],1))'
    coeff=real(fft([value';value[:,N:-1:2]'],1))'
    coefficients=[coeff[:,1]./(2*N) coeff[:,2:N]./N coeff[:,N+1]./(2*N)]
    return coefficients
end


function interpspec2D_FFT(n1::Int64,n2::Int64,fcheby::Array{Float64,2})
coefficients2=zeros(n1,n2)
# First dimension - coefficients for (k)
coefficients2= coeff_fft(fcheby);
# Second dimension - coefficients for (k,l)
coefficients2 = coeff_fft(collect(coefficients2'))';
return coefficients2
end
