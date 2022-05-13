
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
