
##
@inline function getproperty(p::T,f::Symbol) where T<:AbstractNode
    if ~hasfield(T,f)
        data = getfield(p,:data)
        haskey(data,f) ? data[f][getfield(p,:id)] : 0.0
    else
        getfield(p,f)
    end
end

@inline function setproperty!(p::T,f::Symbol,x::Float64) where T<:AbstractNode
    if ~haskey(getfield(p,:data),f)
        n = length(getfield(p,:data)[:w])
        getfield(p,:data)[f] = zeros(n)
    end
    getfield(p,:data)[f][getfield(p,:id)] = x
end

@inline -(n::T,x::NTuple{3,Float64}) where T<:AbstractNode = (n.x-x[1],n.y-x[2],n.z-x[3])
@inline -(x::NTuple{3,Float64},n::T) where T<:AbstractNode = (x[1]-n.x,x[2]-n.y,x[3]-n.z)
@inline +(n::T,x::NTuple{3,Float64}) where T<:AbstractNode = (n.x+x[1],n.y+x[2],n.z+x[3])
@inline +(x::NTuple{3,Float64},n::T) where T<:AbstractNode = (x[1]+n.x,x[2]+n.y,x[3]+n.z)
@inline *(x::NTuple{N,Float64},y::NTuple{N,Float64}) where N = sum(x[i]*y[i] for i in 1:N)

## ----------------- Node ------------------
struct Node<:AbstractNode
    id::Int
    data::Dict{Symbol,Vector{Float64}}
end

## ----------------- SNode ------------------
struct SNode<:AbstractNode
    id::Int
    data::Dict{Symbol,Vector{Float64}}
    index::Vector{Int}
    𝝭::Dict{Symbol,Vector{Float64}}
end

## convert
Node(ξ::SNode) = Node(ξ.id,ξ.data)
SNode(ξ::T,η::SNode) where T<:AbstractNode = SNode(ξ.id,ξ.data,η.index,η.𝝭)
function (η::SNode)(ξ::SNode)
    ξ.index[ξ.id] = η.index[η.id]
    empty!(ξ.𝝭)
    for s in keys(η.𝝭)
        ξ.𝝭[s] = η.𝝭[s]
    end
end
