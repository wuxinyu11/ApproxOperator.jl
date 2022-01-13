
##
@inline getproperty(p::T,f::Symbol) where T<:PhysicalNode = hasfield(T,f) ? getfield(node,f) : getfield(p,:data)[f][getfield(p,:id)]

# ----------------- Node ------------------
struct Node<:AbstractNode
    id::Int
    data::Dict{Symbol,Vector{Float64}}
end
