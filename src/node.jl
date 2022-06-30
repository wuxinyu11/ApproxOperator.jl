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
-(r::RV,a::Float64) = r.v[r.i] - a
-(a::Float64,r::RV) = a - r.v[r.i]
*(r::RV,a::Float64) = r.v[r.i] * a
*(a::Float64,r::RV) = r.v[r.i] * a

getindex(r::RV,i::Int) = r.v[r.i+i]
function setindex!(r::RV,i::Int,x::Float64)
    r.v[r.i+i] = x
end

"""
    GernelNode{S,N}()

GernelNode type
"""
struct GernelNode{S,N}
    index::NamedTuple{S,NTuple{N,Int}}
    data::Dict{Symbol,Tuple{Symbol,Vector{Float64}}}
end

function getproperty(p::GernelNode{S,N},f::Symbol) where {S,N}
    if f∈S
        return getfield(p,:index)[f]
    elseif haskey(getfield(p,:data),f)
        index = getfield(p,:index)
        i,v = getfield(p,:data)[f]
        return RV(index[i],v)
    else
        return 0.0
    end
end

function setproperty!(p::GernelNode,f::Symbol,x::Float64)
    if ~haskey(getfield(p,:data),f)
        i,v = getfield(p,:data)[:x]
        n = length(v)
        getfield(p,:data)[f] = (i,zeros(n))
    end
    index = getfield(p,:index)
    i,v = getfield(p,:data)[f]
    v[index[i]] = x
end
    
"""
    Node
"""
const Node = GernelNode{(:I,),1}
function GernelNode{(:I,),1}(ps::Pair{Symbol,Vector{Float64}}...)
    _,v = ps[1]
    n = length(v)
    data = Dict([s=>(:I,v) for (s,v) in ps])
    return [Node((I=i,),data) for i in 1:n]
end


const SNode = GernelNode{(:g,:G,:s),3}
