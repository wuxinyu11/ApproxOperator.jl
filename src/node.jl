
<<<<<<< HEAD
##
@inline getproperty(p::T,f::Symbol) where T<:PhysicalNode = hasfield(T,f) ? getfield(node,f) : getfield(p,:data)[f][getfield(p,:id)]

# ----------------- Node ------------------
struct Node<:AbstractNode
    id::Int
    data::Dict{Symbol,Vector{Float64}}
end
=======
## PointData
struct Node{N,T}
    i::Int
    data::NTuple{N,Vector{T}}
end
Node(i::Int,datas::AbstractVector...) = Node(i,datas)

getindex(pd::Node,i::Int) = pd.data[i][pd.i]
>>>>>>> d3101db3655d41428b6a1256f1bf9be32a239569
