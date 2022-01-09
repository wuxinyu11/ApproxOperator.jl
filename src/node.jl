@inline +(n::NTuple{3,Float64},m::NTuple{3,Float64}) = (n[1]+m[1], n[2]+m[2], n[3]+m[3])
@inline -(n::NTuple{3,Float64},m::NTuple{3,Float64}) = (n[1]-m[1], n[2]-m[2], n[3]-m[3])
@inline *(c::Float64,n::NTuple{3,Float64}) = (c*n[1], c*n[2], c*n[3])

## PhysicalNode
@inline getproperty(node::T,f::Symbol) where T<:PhysicalNode = hasfield(T,f) ? getfield(node,f) : node.data[f][node.I]

# ----------------- Node ------------------
struct Node<:PhysicalNode
    I::Int
    data::Dict{Symbol,Vector{Float64}}
end

## ParametricNode
@inline getproperty(node::T,f::Symbol) where T<:ParametricNode = hasfield(T,f) ? getfield(node,f) : getdata(node,f)
@inline getdata(node::T,f::Symbol) where T<:ParametricNode = f∈(:ξ,:η,:γ,:w) ? node.data[f][node.G] : node.data[f][node.S]
# ----------------- Node ------------------
struct Point<:ParametricNode
    S::Int
    G::Int
    data::Dict{Symbol,Vector{Float64}}
end
