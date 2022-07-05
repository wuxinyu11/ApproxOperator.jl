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
struct Node<:AbstractNode
    index::Int
    data::Dict{Symbol,Tuple{Int,Vector{Float64}}}
end
const REF = (I=1,)

"""
SNode
"""
struct SNode<:AbstractNode
    index::NTuple{3,Int}
    data::Dict{Symbol,Tuple{Int,Vector{Float64}}}
end
const SREF = (g=1,G=2,s=3)

for (t,ref) in ((:Node,:REF),(:SNode,:SREF))

@eval begin

function getproperty(p::$t,f::Symbol)
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

function setproperty!(p::$t,f::Symbol,x::Float64)
    i,v = getfield(p,:data)[f]
    j = getfield(p,:index)[i]
    v[j] = x
end
#
# function getindex(p::T,f::Symbol)
#     i,v = getfield(p,:data)[f]
#     return RV(i,v)
# end

end
end
