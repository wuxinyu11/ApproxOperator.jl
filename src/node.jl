ParametricCoordinates = Union{Float64,NTuple{2,Float64},NTuple{3,Float64}}
@inline +(n::NTuple{3,Float64},m::NTuple{3,Float64}) = (n[1]+m[1], n[2]+m[2], n[3]+m[3])
@inline -(n::NTuple{3,Float64},m::NTuple{3,Float64}) = (n[1]-m[1], n[2]-m[2], n[3]-m[3])
@inline *(c::Float64,n::NTuple{3,Float64}) = (c*n[1], c*n[2], c*n[3])

@inline getproperty(node::T,f::Symbol) where T<:AbstractNode = getdata(node,Val(f))
@inline function getdata(node::T,::Val{:id}) where T<:AbstractNode
    return getfield(node,:i)
end
@inline function getdata(node::T,::Val{:coordinates}) where T<:AbstractNode
    coordinates = getfield(node,:coordinates)
    i = getfield(node,:i)
    return coordinates[i]
end
@inline function getdata(node::T,::Val{:data}) where T<:AbstractNode
    return getfield(node,:data)
end

## PhysicalNode
@inline function getdata(node::T,::Val{:x}) where T<:PhysicalNode
    coordinates = getfield(node,:coordinates)
    i = getfield(node,:i)
    return coordinates[i][1]
end
@inline function getdata(node::T,::Val{:y}) where T<:PhysicalNode
    coordinates = getfield(node,:coordinates)
    i = getfield(node,:i)
    return coordinates[i][2]
end
@inline function getdata(node::T,::Val{:z}) where T<:PhysicalNode
    coordinates = getfield(node,:coordinates)
    i = getfield(node,:i)
    return coordinates[i][3]
end

for g in (:u, :v, :w, :s₁, :s₂, :s₃)
    @inline function getdata(node::T,::Val{g}) where T<:PhysicalNode
        data = getfield(node,:data)
        i = getfield(node,:i)
        return data[g][i]
    end
end
# ----------------- Node ------------------
struct Node<:PhysicalNode
    i::Int
    coordinates::Vector{NTuple{3,Float64}}
    data::Dict{Symbol,AbstractVector{Float64}}
end

## ParametricNode
@inline function getdata(node::T,::Val{:w}) where T<:ParametricNode
    return getfield(node,:weight)
end
@inline function getdata(node::T,::Val{:ξ}) where T<:ParametricNode
    coordinates = getfield(node,:coordinates)
    j = getfield(node,:j)
    return coordinates[j][1]
end
@inline function getdata(node::T,::Val{:η}) where T<:ParametricNode
    coordinates = getfield(node,:coordinates)
    j = getfield(node,:j)
    return coordinates[j][2]
end
@inline function getdata(node::T,::Val{:γ}) where T<:ParametricNode
    coordinates = getfield(node,:coordinates)
    j = getfield(node,:j)
    return coordinates[j][3]
end

# ----------------- Node ------------------
struct Point<:ParametricNode
    j::Int
    coordinates::Vector{ParametricCoordinates}
    weight::Vector{Float64}
    i::Union{Int,Nothing}
    data::Union{Dict{Symbol,Vector{Float64}},Nothing}
end
