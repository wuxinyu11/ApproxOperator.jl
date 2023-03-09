"""
Returen value (RV)

This is a simple struct
"""
struct RV
    i::Int
    v::Vector{Float64}
end

getindex(r::RV,i::Int) = r.v[r.i+i]
function setindex!(r::RV,x::Float64,i::Int)
    r.v[r.i+i] = x
end


"""
Node
"""
struct Node{T,N}
    index::NTuple{N,Int}
    data::Dict{Symbol,Tuple{Int,Vector{Float64}}}
end

function Base.getproperty(p::Node{T,N},s::Symbol) where {T,N}
    index = getfield(p,:index)
    if haskey(T,s)
        return index[findfirst(x->x==s,keys(T))]
    else
        i,v = getfield(p,:data)[s]
        return v[index[i]]
    end
end
function Base.setproperty!(p::Node,s::Symbol,x::Float64)
    i,v = getfield(p,:data)[s]
    j = getfield(p,:index)[i]
    v[j] = x
end
function Base.getindex(p::Node,f::Symbol)
    i,v = getfield(p,:data)[f]
    j = getfield(p,:index)[i]
    return RV(j,v)
end

+(a::T,b::S) where {T<:Node,S<:Node} = (a.x+b.x,a.y+b.y,a.z+b.z)
-(a::T,b::S) where {T<:Node,S<:Node} = (a.x-b.x,a.y-b.y,a.z-b.z)

function push!(p::Node{T,N},sv::Pair{Symbol,Tuple{Symbol,Vector{Float64}}}) where {T,N}
    (s,v) = sv
    push!(getfield(p,:data),s=>(findfirst(x->x==v[1],keys(T)),v[2]))
end
function push!(p::Node{T,N},sv::Pair{Symbol,Symbol}) where {T,N}
    (s,v) = sv
    push!(getfield(p,:data),s=>(findfirst(x->x==v,keys(T)),zeros(T[v])))
end
push!(ps::Vector{T},sv::Pair) where T<:Node = push!(ps[1],sv)