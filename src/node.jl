
##
@inline getproperty(p::T,f::Symbol) where T<:AbstractNode = hasfield(T,f) ? getfield(p,f) : getfield(p,:data)[f][getfield(p,:id)]
@inline function setproperty!(p::T,f::Symbol,x::Float64) where T<:AbstractNode
    getfield(p,:data)[f][getfield(p,:id)] = x
end

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
