
# ##
# @inline function getproperty(p::T,f::Symbol) where T<:AbstractNode
#     if ~hasfield(T,f)
#         data = getfield(p,:data)
#         haskey(data,f) ? data[f][getfield(p,:id)] : 0.0
#     else
#         getfield(p,f)
#     end
# end

# @inline function setproperty!(p::T,f::Symbol,x::Float64) where T<:AbstractNode
#     if ~haskey(getfield(p,:data),f)
#         n = length(getfield(p,:data)[:w])
#         getfield(p,:data)[f] = zeros(n)
#     end
#     getfield(p,:data)[f][getfield(p,:id)] = x
# end

# @inline -(n::T,x::NTuple{3,Float64}) where T<:AbstractNode = (n.x-x[1],n.y-x[2],n.z-x[3])
# @inline -(x::NTuple{3,Float64},n::T) where T<:AbstractNode = (x[1]-n.x,x[2]-n.y,x[3]-n.z)
# @inline +(n::T,x::NTuple{3,Float64}) where T<:AbstractNode = (n.x+x[1],n.y+x[2],n.z+x[3])
# @inline +(x::NTuple{3,Float64},n::T) where T<:AbstractNode = (x[1]+n.x,x[2]+n.y,x[3]+n.z)
# @inline *(x::NTuple{N,Float64},y::NTuple{N,Float64}) where N = sum(x[i]*y[i] for i in 1:N)
# @inline /(x::NTuple{3,Float64},c::Float64) = (x[1]/c,x[2]/c,x[3]/c)

# ## ----------------- Node ------------------
# struct Node<:AbstractNode
#     id::Int
#     data::Dict{Symbol,Vector{Float64}}
# end

struct Node{N}
    index::NTuple{N,Int}
    data::Dict{Symbol,Tuple{Int,Vector{Float64}}}
end

function getproperty(p::Node,f::Symbol)
    if haskey(getfield(p,:data),f)
        index = getfield(p,:index)
        i,v = getfield(p,:data)[f]
        return v[index[i]]
    else
        return 0.0
    end
end

function setproperty!(p::Node,f::Symbol,x::Float64)
    if ~haskey(getfield(p,:data),f)
        i,v = getfield(p,:data)[:x]
        n = length(v)
        getfield(p,:data)[f] = (i,zeros(n))
    end
    index = getfield(p,:index)
    i,v = getfield(p,:data)[f]
    v[index[i]] = x
end

function getindex(p::Node,f::Symbol,j::Int)
    index = getfield(p,:index)
    i,v = getfield(p,:data)[f]
    return v[index[i]+j]
end

function setindex!(p::Node,x::Float64,f::Symbol,j::Int)
    index = getfield(p,:index)
    i,v = getfield(p,:data)[f]
    v[index[i]+j] = x
end