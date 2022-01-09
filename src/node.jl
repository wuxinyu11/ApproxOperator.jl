@inline +(n::NTuple{3,Float64},m::NTuple{3,Float64}) = (n[1]+m[1], n[2]+m[2], n[3]+m[3])
@inline -(n::NTuple{3,Float64},m::NTuple{3,Float64}) = (n[1]-m[1], n[2]-m[2], n[3]-m[3])
@inline *(c::Float64,n::NTuple{3,Float64}) = (c*n[1], c*n[2], c*n[3])

@inline getproperty(node::T,f::Symbol) where T<:AbstractNode = hasfield(T,f) ? getfield(node,f) : getdata(node,f)

## PhysicalNode
# @inline getdata(node::T,f::Symbol) where T<:PhysicalNode = f===:coordinates ? (node.x,node.y,node.z) : node.data[f][node.I]
@inline getdata(node::T,f::Symbol) where T<:PhysicalNode = f==:coordinates ? getdata(node,Val(f)) : node.data[f][node.I]
@inline getdata(node::T,f::Val{:coordinates}) where T<:PhysicalNode = (node.x,node.y,node.z)

# ----------------- Node ------------------
struct Node<:PhysicalNode
    I::Int
    data::Dict{Symbol,Vector{Float64}}
end

## ParametricNode
@inline getdata(node::T,f::Symbol) where T<:ParametricNode = f∈(:coordinates,:ξ,:η,:γ,:w) ? (f===:coordinates ? getdata(node,Val(f)) : node.paradata[f][node.G]) : node.physdata[f][node.S]
@inline getdata(node::T,::Val{:coordinates}) where T<:ParametricNode = haskey(node.paradata,:γ) ? (node.ξ,node.η,node.γ) : (haskey(node.paradata,:η) ? (node.ξ,node.η) : node.ξ)

# ----------------- Node ------------------
struct Point<:ParametricNode
    G::Int
    paradata::Union{Dict{Symbol,Vector{Float64}},Nothing}
    S::Union{Int,Nothing}
    physdata::Union{Dict{Symbol,Vector{Float64}},Nothing}
end
