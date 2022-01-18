
##
@inline getproperty(p::T,f::Symbol) where T<:AbstractNode = hasfield(T,f) ? getfield(p,f) : getfield(p,:data)[f][getfield(p,:id)]
@inline function setproperty!(p::T,f::Symbol,x::Float64) where T<:AbstractNode
    getfield(p,:data)[f][getfield(p,:id)] = x
end

@inline -(n::T,x::NTuple{3,Float64}) where T<:AbstractNode = (n.x-x[1],n.y-x[2],n.z-x[3])
@inline -(x::NTuple{3,Float64},n::T) where T<:AbstractNode = (x[1]-n.x,x[2]-n.y,x[3]-n.z)
@inline +(n::T,x::NTuple{3,Float64}) where T<:AbstractNode = (n.x+x[1],n.y+x[2],n.z+x[3])
@inline +(x::NTuple{3,Float64},n::T) where T<:AbstractNode = (x[1]+n.x,x[2]+n.y,x[3]+n.z)

# ----------------- Node ------------------
struct Node<:AbstractNode
    id::Int
    data::Dict{Symbol,Vector{Float64}}
end

## Meshfree module
# ----------------- MFNode ------------------
struct MFNode{𝒑,𝑠,𝜙}<:AbstractNode
    id::Int
    data::Dict{Symbol,Vector{Float64}}
    type::Tuple{Val{𝒑},Val{𝑠},Val{𝜙}}
end
