## Actions for general functions
# @inline get_shape_functions(ap::Approximator,ξ::ParametricCoordinates,gs::Val...) = (get_shape_functions(ap,ξ,g) for g in gs)
setindex!(ap::Approximator,i::Int) = ap.𝓒[i]
length(ap::Approximator) = length(ap.𝓒)

## AbstractPoi
get_jacobe(::AbstractPoi,::Float64) = 1.
get_coordinates(ap::AbstractPoi,::Float64) = ap.𝓧[ap.𝓒[1]].coordinates

# --------------- Poi1 ---------------
struct Poi1 <: AbstractPoi
    𝓒::Int
    𝓧::Vector{PhysicalNode}
    𝓖::Vector{ParametricNode}
end

# constructions of Poi1
Poi1(𝓒::Int,𝓧::Vector{T};𝓖::Symbol=:PoiGI1) where T<:PhysicalNode = Poi1(𝓒,𝓧,QuadratureRule[𝓖])

# actions of Poi1
get_shape_functions(::Poi1,::Float64,::Val{:∂1}) = 1.


## AbstractSeg
get_jacobe(ap::AbstractSeg,::Float64) = 0.5*ap.L
function get_coordinates(ap::AbstractSeg,ξ::Float64)
    N₁ = (1.0-ξ)*0.5
    N₂ = (1.0+ξ)*0.5
    return N₁*ap.𝓧[ap.𝓒[1]].coordinates + N₂*ap.𝓧[ap.𝓒[2]].coordinates
end
function get_coordinates(ap1::AbstractSeg,ap2::AbstractPoi,::Float64)
    c₁ = findfirst(x -> x == ap2.𝓒[1], ap1.𝓒)
    return (c₁ == 1 ? -1. : 1.)
end
function get_normal(ap::AbstractSeg)
    x₁ = ap.𝓧[ap.𝓒[1]].coordinates[1]
    y₁ = ap.𝓧[ap.𝓒[1]].coordinates[2]
    x₂ = ap.𝓧[ap.𝓒[2]].coordinates[1]
    y₂ = ap.𝓧[ap.𝓒[2]].coordinates[2]
    L = ap.L
    return (y₂-y₁)/L,(x₁-x₂)/L,0.
end
function get_normal(ap1::AbstractSeg,ap2::AbstractPoi)
    c₁ = findfirst(x -> x == ap2.𝓒[1], ap1.𝓒)
    return (c₁ == 1 ? (-1.,0.,0.) : (1.,0.,0.))
end

# --------------- Seg2 ---------------
mutable struct Seg2 <: AbstractSeg
    𝓒::Vector{Int}
    𝓧::Vector{PhysicalNode}
    𝓖::Vector{ParametricNode}
    L::Float64
end

# constructions of Seg2
function Seg2(𝓒::Vector{Int},𝓧::Vector{T};𝓖::Symbol=:SegGI2) where T<:PhysicalNode
    x₁ = 𝓧[𝓒[1]].coordinates[1]
    y₁ = 𝓧[𝓒[1]].coordinates[2]
    x₂ = 𝓧[𝓒[2]].coordinates[1]
    y₂ = 𝓧[𝓒[2]].coordinates[2]
    L = ((x₂-x₁)^2+(y₂-y₁)^2)^0.5
    𝓖 = QuadratureRule[𝓖]
    return Seg2(𝓒,𝓧,𝓖,L)
end

# actions of Seg2
get_shape_functions(::Seg2,ξ::Float64,::Val{:∂1}) = ((1.0-ξ)*0.5,(1.0+ξ)*0.5)
function get_shape_functions(ap::Seg2,::Float64,::Val{:∂x})
    x₁ = ap.𝓧[ap.𝓒[1]].coordinates[1]
    x₂ = ap.𝓧[ap.𝓒[2]].coordinates[1]
    return (-1.0/(x₂-x₁),1.0/(x₂-x₁))
end
get_shape_functions(::Seg2,ξ::Float64,::Val{:∂y}) = (0.,0.)
get_shape_functions(::Seg2,ξ::Float64,::Val{:∂z}) = (0.,0.)

## AbstractTri
get_jacobe(ap::AbstractTri,ξ::NTuple{2,Float64}) = ap.A
function get_coordinates(ap::AbstractTri,ξ::NTuple{2,Float64})
    N₁ = ξ[1]
    N₂ = ξ[2]
    N₃ = 1-ξ[1]-ξ[2]
    return N₁*ap.𝓧[ap.𝓒[1]].coordinates + N₂*ap.𝓧[ap.𝓒[2]].coordinates + N₃*ap.𝓧[ap.𝓒[3]].coordinates
end
function get_coordinates(ap1::AbstractTri,ap2::AbstractSeg,ξ::NTuple{2,Float64})
    c₁ = findfirst(x -> x == ap2.𝓒[1], ap1.𝓒)
    c₂ = findfirst(x -> x == ap2.𝓒[2], ap1.𝓒)
    return SVector{2,Float64}((1-ξ[1])/2*(c₁ == 1) + (1+ξ[1])/2*(c₂ == 1),
                              (1-ξ[1])/2*(c₁ == 2) + (1+ξ[1])/2*(c₂ == 2))
end
function get_normal(ap::AbstractTri)
    x₁ = ap.𝓧[ap.𝓒[1]].coordinates[1]
    y₁ = ap.𝓧[ap.𝓒[1]].coordinates[2]
    z₁ = ap.𝓧[ap.𝓒[1]].coordinates[3]
    x₂ = ap.𝓧[ap.𝓒[2]].coordinates[1]
    y₂ = ap.𝓧[ap.𝓒[2]].coordinates[2]
    z₂ = ap.𝓧[ap.𝓒[2]].coordinates[3]
    x₃ = ap.𝓧[ap.𝓒[3]].coordinates[1]
    y₃ = ap.𝓧[ap.𝓒[3]].coordinates[2]
    z₃ = ap.𝓧[ap.𝓒[3]].coordinates[3]
    A₁ = 0.5*(y₁*z₂+y₂*z₃+y₃*z₁-y₂*z₁-y₃x*z₂-y₁*z₃)
    A₂ = 0.5*(z₁*x₂+z₂*x₃+z₃*x₁-z₂*x₁-z₃x*x₂-z₁*x₃)
    A₃ = 0.5*(x₁*y₂+x₂*y₃+x₃*y₁-x₂*y₁-x₃x*y₂-x₁*y₃)
    A = ap.A
    return Ax/A,Ay/A,Az/A
end
function get_normal(ap1::AbstractTri,ap2::AbstractSeg)
    c₁ = findfirst(x -> x == ap2.𝓒[1], ap1.𝓒)
    c₂ = findfirst(x -> x == ap2.𝓒[2], ap1.𝓒)
    x₁ = ap1.𝓧[ap1.𝓒[c₁]].coordinates[1]
    y₁ = ap1.𝓧[ap1.𝓒[c₁]].coordinates[2]
    x₂ = ap1.𝓧[ap1.𝓒[c₂]].coordinates[1]
    y₂ = ap1.𝓧[ap1.𝓒[c₂]].coordinates[2]
    L = ap2.L
    return (y₂-y₁)/L,(x₁-x₂)/L,0.
end

# --------------- Tri3 ---------------
# Constant strain triangular Approximator (CST)
struct Tri3 <: AbstractTri
    𝓒::Vector{Int}
    𝓧::Vector{PhysicalNode}
    𝓖::Vector{ParametricNode}
    A::Float64
end

# constructions
function Tri3(𝓒::Vector{Int},𝓧::Vector{T};𝓖::Symbol=:TriGI3) where T<:PhysicalNode
    x₁ = 𝓧[𝓒[1]].coordinates[1]
    y₁ = 𝓧[𝓒[1]].coordinates[2]
    z₁ = 𝓧[𝓒[1]].coordinates[3]
    x₂ = 𝓧[𝓒[2]].coordinates[1]
    y₂ = 𝓧[𝓒[2]].coordinates[2]
    z₂ = 𝓧[𝓒[2]].coordinates[3]
    x₃ = 𝓧[𝓒[3]].coordinates[1]
    y₃ = 𝓧[𝓒[3]].coordinates[2]
    z₃ = 𝓧[𝓒[3]].coordinates[3]
    A₁ = 0.5*(y₁*z₂+y₂*z₃+y₃*z₁-y₂*z₁-y₃*z₂-y₁*z₃)
    A₂ = 0.5*(z₁*x₂+z₂*x₃+z₃*x₁-z₂*x₁-z₃*x₂-z₁*x₃)
    A₃ = 0.5*(x₁*y₂+x₂*y₃+x₃*y₁-x₂*y₁-x₃*y₂-x₁*y₃)
    A = (A₁^2 + A₂^2 + A₃^2)^0.5
    𝓖 = QuadratureRule[𝓖]
    return Tri3(𝓒,𝓧,𝓖,A)
end

# actions
get_shape_functions(ap::Tri3,ξ::NTuple{2,Float64},::Val{:∂1}) = SVector{3,Float64}(ξ[1],ξ[2],1-ξ[1]-ξ[2])
function get_shape_functions(ap::Tri3,ξ::NTuple{2,Float64},::Val{:∂x})
    y₁ = ap.𝓧[ap.𝓒[1]].coordinates[2]
    y₂ = ap.𝓧[ap.𝓒[2]].coordinates[2]
    y₃ = ap.𝓧[ap.𝓒[3]].coordinates[2]
    A = ap.A
    return ((y₂-y₃)/(2^A),(y₃-y₁)/(2*A),(y₁-y₃)/(2*A))
end
function get_shape_functions(ap::Tri3,ξ::NTuple{2,Float64},::Val{:∂y})
    x₁ = ap.𝓧[ap.𝓒[1]].coordinates[1]
    x₂ = ap.𝓧[ap.𝓒[2]].coordinates[1]
    x₃ = ap.𝓧[ap.𝓒[3]].coordinates[1]
    A = ap.A
    return ((x₃-x₂)/(2*A),(x₁-x₃)/(2*A),(x₂-x₁)/(2*A))
end
get_shape_functions(ap::Tri3,ξ::NTuple{2,Float64},::Val{:∂z}) = (0.,0.,0.)

## AbstractQuad
function get_jacobe(ap::AbstractQuad,ξ::NTuple{2,Float64})
    J₁₁,J₂₁,J₁₂,J₂₂ = get_jacobe_matrix(ap,ξ)
    return J₁₁*J₂₂-J₂₁*J₁₂
end
function get_jacobe_matrix(ap::AbstractQuad,ξ::NTuple{2,Float64})
    x₁ = ap.𝓧[ap.𝓒[1]].coordinates[1]
    x₂ = ap.𝓧[ap.𝓒[2]].coordinates[1]
    x₃ = ap.𝓧[ap.𝓒[3]].coordinates[1]
    x₄ = ap.𝓧[ap.𝓒[4]].coordinates[1]
    y₁ = ap.𝓧[ap.𝓒[1]].coordinates[2]
    y₂ = ap.𝓧[ap.𝓒[2]].coordinates[2]
    y₃ = ap.𝓧[ap.𝓒[3]].coordinates[2]
    y₄ = ap.𝓧[ap.𝓒[4]].coordinates[2]
    ∂N₁∂ξ,∂N₂∂ξ,∂N₃∂ξ,∂N₄∂ξ = get_shape_functions(ap,ξ,Val(:∂ξ))
    ∂N₁∂η,∂N₂∂η,∂N₃∂η,∂N₄∂η = get_shape_functions(ap,ξ,Val(:∂η))
    J₁₁ = ∂N₁∂ξ*x₁ + ∂N₂∂ξ*x₂ + ∂N₃∂ξ*x₃ + ∂N₄∂ξ*x₄
    J₁₂ = ∂N₁∂η*x₁ + ∂N₂∂η*x₂ + ∂N₃∂η*x₃ + ∂N₄∂η*x₄
    J₂₁ = ∂N₁∂ξ*y₁ + ∂N₂∂ξ*y₂ + ∂N₃∂ξ*y₃ + ∂N₄∂ξ*y₄
    J₂₂ = ∂N₁∂η*y₁ + ∂N₂∂η*y₂ + ∂N₃∂η*y₃ + ∂N₄∂η*y₄
    return J₁₁,J₂₁,J₁₂,J₂₂
end
function get_coordinates(ap::AbstractQuad,ξ::NTuple{2,Float64})
    N₁,N₂,N₃,N₄ = get_shape_functions(ap,ξ,Val(:∂1))
    return N₁*ap.𝓧[ap.𝓒[1]].coordinates + N₂*ap.𝓧[ap.𝓒[2]].coordinates + N₃*ap.𝓧[ap.𝓒[3]].coordinates + N₄*ap.𝓧[ap.𝓒[4]].coordinates
end
function get_coordinates(ap1::AbstractQuad,ap2::AbstractSeg,ξ::NTuple{2,Float64})
    c₁ = findfirst(x -> x == ap2.𝓒[1], ap1.𝓒)
    c₂ = findfirst(x -> x == ap2.𝓒[2], ap1.𝓒)
    return SVector{2,Float64}((1-ξ[1])/2*(c₁ == 1)*(c₂ == 2) + (1+ξ[1])/2*(c₁ == 3)*(c₂ == 4),
                              (1-ξ[1])/2*(c₁ == 2)*(c₂ == 3) + (1+ξ[1])/2*(c₁ == 4)*(c₂ == 1))
end

# --------------- Quad ---------------
mutable struct Quad <: AbstractQuad
    𝓒::Vector{Int}
    𝓧::Vector{PhysicalNode}
    𝓖::Vector{ParametricNode}
end
# constructions
function Quad(𝓒::Vector{Int},𝓧::Vector{T};𝓖::Symbol=:QuadGI2) where T<:PhysicalNode
    𝓖 = QuadratureRule[𝓖]
    return Quad(𝓒,𝓧,𝓖)
end

# actions
function get_shape_functions(ap::Quad,ξ::NTuple{2,Float64},::Val{:∂1})
    N₁ = 0.25*(1-ξ[1])*(1-ξ[2])
    N₂ = 0.25*(1+ξ[1])*(1-ξ[2])
    N₃ = 0.25*(1+ξ[1])*(1+ξ[2])
    N₄ = 0.25*(1-ξ[1])*(1+ξ[2])
    return (N₁,N₂,N₃,N₄)
end
function get_shape_functions(ap::Quad,ξ::NTuple{2,Float64},::Val{:∂ξ})
    ∂N₁∂ξ = - 0.25*(1-ξ[2])
    ∂N₂∂ξ =   0.25*(1-ξ[2])
    ∂N₃∂ξ =   0.25*(1+ξ[2])
    ∂N₄∂ξ = - 0.25*(1+ξ[2])
    return (∂N₁∂ξ,∂N₂∂ξ,∂N₃∂ξ,∂N₄∂ξ)
end
function get_shape_functions(ap::Quad,ξ::NTuple{2,Float64},::Val{:∂η})
    ∂N₁∂η = - 0.25*(1-ξ[1])
    ∂N₂∂η = - 0.25*(1+ξ[1])
    ∂N₃∂η =   0.25*(1+ξ[1])
    ∂N₄∂η =   0.25*(1-ξ[1])
    return (∂N₁∂η,∂N₂∂η,∂N₃∂η,∂N₄∂η)
end
function get_shape_functions(ap::Quad,ξ::NTuple{2,Float64},::Val{:∂x},::Val{:∂y})
    x₁ = ap.𝓧[ap.𝓒[1]].coordinates[1]
    x₂ = ap.𝓧[ap.𝓒[2]].coordinates[1]
    x₃ = ap.𝓧[ap.𝓒[3]].coordinates[1]
    x₄ = ap.𝓧[ap.𝓒[4]].coordinates[1]
    y₁ = ap.𝓧[ap.𝓒[1]].coordinates[2]
    y₂ = ap.𝓧[ap.𝓒[2]].coordinates[2]
    y₃ = ap.𝓧[ap.𝓒[3]].coordinates[2]
    y₄ = ap.𝓧[ap.𝓒[4]].coordinates[2]
    ∂N₁∂ξ,∂N₂∂ξ,∂N₃∂ξ,∂N₄∂ξ = get_shape_functions(ap,ξ,Val(:∂ξ))
    ∂N₁∂η,∂N₂∂η,∂N₃∂η,∂N₄∂η = get_shape_functions(ap,ξ,Val(:∂η))
    ∂x∂ξ = ∂N₁∂ξ*x₁ + ∂N₂∂ξ*x₂ + ∂N₃∂ξ*x₃ + ∂N₄∂ξ*x₄
    ∂x∂η = ∂N₁∂η*x₁ + ∂N₂∂η*x₂ + ∂N₃∂η*x₃ + ∂N₄∂η*x₄
    ∂y∂ξ = ∂N₁∂ξ*y₁ + ∂N₂∂ξ*y₂ + ∂N₃∂ξ*y₃ + ∂N₄∂ξ*y₄
    ∂y∂η = ∂N₁∂η*y₁ + ∂N₂∂η*y₂ + ∂N₃∂η*y₃ + ∂N₄∂η*y₄
    detJ = ∂x∂ξ*∂y∂η - ∂x∂η*∂y∂ξ
    ∂ξ∂x =   ∂y∂η/detJ
    ∂η∂x = - ∂y∂ξ/detJ
    ∂ξ∂y = - ∂x∂η/detJ
    ∂η∂y =   ∂x∂ξ/detJ
    ∂N₁∂x = ∂N₁∂ξ*∂ξ∂x + ∂N₁∂η*∂η∂x
    ∂N₂∂x = ∂N₂∂ξ*∂ξ∂x + ∂N₂∂η*∂η∂x
    ∂N₃∂x = ∂N₃∂ξ*∂ξ∂x + ∂N₃∂η*∂η∂x
    ∂N₄∂x = ∂N₄∂ξ*∂ξ∂x + ∂N₄∂η*∂η∂x
    ∂N₁∂y = ∂N₁∂ξ*∂ξ∂y + ∂N₁∂η*∂η∂y
    ∂N₂∂y = ∂N₂∂ξ*∂ξ∂y + ∂N₂∂η*∂η∂y
    ∂N₃∂y = ∂N₃∂ξ*∂ξ∂y + ∂N₃∂η*∂η∂y
    ∂N₄∂y = ∂N₄∂ξ*∂ξ∂y + ∂N₄∂η*∂η∂y
    return (∂N₁∂x,∂N₂∂x,∂N₃∂x,∂N₄∂x),(∂N₁∂y,∂N₂∂y,∂N₃∂y,∂N₄∂y)
end

get_shape_functions(ap::Quad,ξ::NTuple{2,Float64},::Val{:∂z}) = (0.,0.,0.,0.)
@inline get_shape_functions(ap::Quad,ξ::NTuple{2,Float64},::Val{:∂1},::Val{:∂x},::Val{:∂y}) = get_shape_functions(ap,ξ,Val(:∂1)),get_shape_functions(ap,ξ,Val(:∂x),Val(:∂y))...
@inline get_shape_functions(ap::Quad,ξ::NTuple{2,Float64},::Val{:∂1},::Val{:∂x},::Val{:∂y},::Val{:∂z}) = get_shape_functions(ap,ξ,Val(:∂1)),get_shape_functions(ap,ξ,Val(:∂x),Val(:∂y))...,get_shape_functions(ap,ξ,Val(:∂z))
