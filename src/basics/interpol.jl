#License is MIT: https://github.com/sbadred/LTEI_TA.jl/blob/99b988ec2d84266e51a5a9b6a5acaf190c26e019/LICENSE

"""Get the optimal number of interpolation points for a given
smooth functions
##############################################################
Input:
-cluster: Computational box
-t : quadrature nodes
-maxi: maximum number of interpolation points you need

Output:
-List of interpolation nodes
/***************************************************************************************
*    Title: chebfun
*    Date: 2017
#    Version: CHEBFUN v5.7.0 RELEASE
*    Availability: https://github.com/chebfun/chebfun/archive/master.zip
*
***************************************************************************************/

##############################################################
"""


using MATLAB
mat"addpath('/Users/sbadredd/Desktop/low-rank-two-electron-integrals/SOLVER/chebfun-chebfun-c07b658')"


function interpolation(max_N1::Int64,n1::Array{Int64,1};flag::Bool=true)
    idx=findall(x->x>max_N1, n1);
    if !isempty(idx)
        n1[idx] .=max_N1;
    else
        if flag
            n1=max_N1 .*Int.(ones(length(n1)))
        else
            return n1
        end
    end

    return n1
end

function interpol(cluster ,t::Array{Float64,1},maxi::Int64;flag::Bool=true)
    n1=zeros(Float64,length(t));
    try
        @inbounds for i in 1:length(t)
            boite=convert.(AbstractFloat,[cluster cluster])
            mat"f=chebfun2(@(x,y) exp(-$t($i)^2.*(x-y).^2),$boite);n=rank(f);"
            n1[i]=@mget n
            println(i)
        end
        n1=interpolation(maxi,Int.(n1),flag=flag)
    catch
        println("aborted")
    end
    return n1

end

function number_INT(boite::AbstractArray,
    t::Array{Float64,1},maxi::Int64;Timelimit::Float64=5,flag::Bool=true)
    task = @async(Int.(interpol(boite,t,maxi,flag=flag)))  # run everything via myfunc()
    # run everything via myfunc()
    println("task is started: ", istaskstarted(task))
    println("Task is done: ", istaskdone(task))
    sleep(Timelimit) # wait a while ...
    #@show "i slept"
    if !istaskdone(task)
        println("Killing task.")
        @async(Base.throwto(task, DivideError()));
        n1=max_N1 .*Int.(ones(length(t)))
        return n1
    end
    return fetch(task)

end

function interpol(cluster::Array{clus,1},t::Array{Float64,1},maxi::Int64;flag::Bool=true)
    n1=zeros(Float64,length(t),length(cluster));
    try
        @inbounds for j in 1:length(cluster)
            @inbounds for i in 1:length(t)
                boite=convert.(AbstractFloat,[cluster[j].boite cluster[j].boite])
                mat"f=chebfun2(@(x,y) exp(-$t($i)^2.*(x-y).^2),$boite);n=rank(f);"
                #f=Fun((x,y) ->exp(-t[i]^2 .*(x-y).^2),Interval(cluster[j].boite[1],cluster[j].boite[2])^2)
                #n1[i,j]=ncoefficients(f);
                n1[i,j]=@mget n
            end
        end
        n1=interpolation(maxi,Int.(n1),flag=flag)
    catch
        println("aborted")
    end
    return n1
end
