
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

# ----------------- MFNode ------------------
struct MFNode{𝒑,𝑠,𝜙}<:AbstractNode
    id::Int
    data::Dict{Symbol,Vector{Float64}}
    type::Tuple{Val{𝒑},Val{𝑠},Val{𝜙}}
end

## Quadrature Point
# push!(ap::A,s::Symbol) where A<:AbstractSeg = push!(ap,Val(s))
# function push!(ap::A,data::Dict{Symbol,Vector{Float64}},v::NTuple{2,Float64}) where A<:AbstractSeg
#     ξ = data[:ξ]
#     w = data[:w]
#     push!(ξ,v[1])
#     push!(w,v[2])
#     id = length(w)
#     push!(ap.𝓖,Node(id,data))
# end
#
# function push!(aps::Vector{A},data)
