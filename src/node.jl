"""
Returen value (RV)

This is a simple struct
"""
struct RV
    i::Int
    v::Vector{Float64}
end

"""
+, -, * of type RV
"""
+(r::RV,a::Float64) = r.v[r.i] + a
+(a::Float64,r::RV) = r.v[r.i] + a
+(r::RV,s::RV) = r.v[r.i] + s.v[s.i]
-(r::RV,a::Float64) = r.v[r.i] - a
-(a::Float64,r::RV) = a - r.v[r.i]
-(r::RV,s::RV) = r.v[r.i] - s.v[s.i]
*(r::RV,a::Float64) = r.v[r.i] * a
*(a::Float64,r::RV) = r.v[r.i] * a
*(r::RV,s::RV) = r.v[r.i] * s.v[s.i]

getindex(r::RV,i::Int) = r.v[r.i+i]
function setindex!(r::RV,x::Float64,i::Int)
    r.v[r.i+i] = x
end

"""
Node
"""
+(a::T,b::S) where {T<:AbstractNode,S<:AbstractNode} = (a.x+b.x,a.y+b.y,a.z+b.z)
-(a::T,b::S) where {T<:AbstractNode,S<:AbstractNode} = (a.x-b.x,a.y-b.y,a.z-b.z)

struct Node<:AbstractNode
    index::Int
    data::Dict{Symbol,Tuple{Int,Vector{Float64}}}
end
const REF = (𝐼=1,)

function Node(ss::Pair{Symbol,Vector{Float64}}...)
    data = Dict([s=>(1,v) for (s,v) in ss])
    _,v = ss[1]
    n = length(v)
    return [Node(i,data) for i in 1:n]
end

"""
SNode
"""
struct SNode<:AbstractNode
    index::NTuple{4,Int}
    data::Dict{Symbol,Tuple{Int,Vector{Float64}}}
end
const SREF = (𝑔=1,𝐺=2,𝐶=3,𝑠=4)

"""
GNode
"""
struct GNode<:AbstractNode
    index::NTuple{2,Int}
    data::Dict{Symbol,Tuple{Int,Vector{Float64}}}
end
const GREF = (𝑖=1,𝐼=2)

for (t,ref) in ((:Node,:REF),(:SNode,:SREF),(:GNode,:GREF))
    @eval begin
        function $t(data::Dict{Symbol,Tuple{Int,Vector{Float64}}},𝐼s::Pair{Symbol,Int}...)
            index = (haskey(𝐼s,s) ? 𝐼s[s] : 0 for (s,i) in $ref)
            return $t(index,data)
        end

        function Base.getproperty(p::$t,f::Symbol)
            if haskey($ref,f)
                return getfield(p,:index)[$ref[f]]
            elseif haskey(getfield(p,:data),f)
                i,v = getfield(p,:data)[f]
                j = getfield(p,:index)[i]
                return v[j]
            else
                return 0.0
            end
        end

        function Base.setproperty!(p::$t,f::Symbol,x::Float64)
            i,v = getfield(p,:data)[f]
            j = getfield(p,:index)[i]
            v[j] = x
        end

        function Base.getindex(p::$t,f::Symbol)
            i,v = getfield(p,:data)[f]
            j = getfield(p,:index)[i]
            return RV(j,v)
        end

    end
end

"""
push!(node<:AbstractNode)
"""
function push!(n::AbstractNode,svs::Pair{Symbol,Vector{Float64}}...;index::Int=1)
    for (s,v) in svs
        push!(getfield(n,:data),s=>(index,v))
    end
end
function push!(ns::Vector{N},svs::Pair{Symbol,Vector{Float64}}...;index::Int=1) where N<:AbstractNode
    push!(ns[1],svs...,index=index)
end
