```
Node
```
struct Node{S,N}
    index::NamedTuple{S,NTuple{N,Int}}
    data::Dict{Symbol,Tuple{Symbol,Vector{Float64}}}
end

function getproperty(p::Node{S,N},f::Symbol) where {S,N}
    if f∈S
        return getfield(p,:index)[f]
    elseif haskey(getfield(p,:data),f)
        index = getfield(p,:index)
        i,v = getfield(p,:data)[f]
        return v[index[i]]
    else
        return 0.0
    end
end

function setproperty!(p::Node,f::Symbol,x::Float64)
    if ~haskey(getfield(p,:data),f)
        i,v = getfield(p,:data)[:x]
        n = length(v)
        getfield(p,:data)[f] = (i,zeros(n))
    end
    index = getfield(p,:index)
    i,v = getfield(p,:data)[f]
    v[index[i]] = x
end

function getindex(p::Node,f::Symbol,j::Int)
    index = getfield(p,:index)
    i,v = getfield(p,:data)[f]
    return v[index[i]+j]
end

function setindex!(p::Node,x::Float64,f::Symbol,j::Int)
    index = getfield(p,:index)
    i,v = getfield(p,:data)[f]
    v[index[i]+j] = x
end

```
Special Node
```
FiniteNode = Node{(:I),1}
QuadratureNode = Node{(:g,:G,:s),3}