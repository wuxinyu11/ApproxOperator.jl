@inline +(n::NTuple{3,Float64},m::NTuple{3,Float64}) = (n[1]+m[1], n[2]+m[2], n[3]+m[3])
@inline -(n::NTuple{3,Float64},m::NTuple{3,Float64}) = (n[1]-m[1], n[2]-m[2], n[3]-m[3])
@inline *(c::Float64,n::NTuple{3,Float64}) = (c*n[1], c*n[2], c*n[3])

##
struct Node<:AbstractNode
    coordinates::Coordinates
end
@inline getproperty(node::Node,f::Symbol) = getdata(node,Val(f))
@inline function getproperty(node::Node,::Val{:x})
    coordinates = getfield(node,:coordinates)
    return coordinates.x
end
@inline function getproperty(node::Node,::Val{:y})
    coordinates = getfield(node,:coordinates)
    return coordinates.y
end
@inline function getproperty(node::Node,::Val{:z})
    coordinates = getfield(node,:coordinates)
    return coordinates.z
end

struct Point<:AbstractNode
    cooridnates::ParametricCoordinates
end
