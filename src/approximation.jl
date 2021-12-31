## Actions for general functions
@inline get_global_indice(ap::Approximator,i::Int) = ap.id[i]
@inline get_number_of_indices(ap::Approximator) = length(ap.id)
@inline get_local_node(ap::Approximator,i::Int) = ap.nd[ap.id[i]].x
@inline get_shape_functions(ap::Approximator,ξ::AbstractVector{Float64},gs::Val...) = (get_shape_functions(ap,ξ,g) for g in gs)

## AbstractPoi
get_jacobe(::AbstractPoi,::AbstractVector{Float64}) = 1.
get_coordinates(ap::AbstractPoi,::AbstractVector{Float64}) = 1.0*ap.nd[ap.id[1]]

# --------------- Poi1 ---------------
struct Poi1 <: AbstractPoi
    𝒙::Vector{PhysicalNode}
    𝓘::Int
    𝓖::Vector{ParametricNode}
end

# constructions of Poi1
Poi1(nd::Vector{Node},d::Int;qw::Symbol=:PoiGI1) = Poi1(nd,id,QuadratureRule[qw])

# actions of Poi1
get_shape_functions(::Poi1,::AbstractVector{Float64},::Val{:∂1}) = 1.
get_shape_functions(::Poi1,::AbstractVector{Float64},::Val{:∂x}) = 1.


## AbstractSeg
get_jacobe(ap::AbstractSeg,::AbstractVector{Float64}) = ap.norm/2
function get_coordinates(ap::AbstractSeg,ξ::AbstractVector{Float64})
    N1 = (1.0-ξ[1])*0.5
    N2 = (1.0+ξ[1])*0.5
    return N1*ap.nd[ap.id[1]].x + N2*ap.nd[ap.id[2]].x
end
function get_coordinates(ap1::AbstractSeg,ap2::AbstractPoi,::AbstractVector{Float64})
    id₁ = findfirst(x -> x == ap2.id[1], ap1.id)
    return (id₁ == 1 ? -1. : 1.)
end
function get_normal(ap::AbstractSeg)
    L = ap.nm
    x1 = ap.nd[ap.id[1]].x[1]
    y1 = ap.nd[ap.id[1]].x[2]
    x2 = ap.nd[ap.id[2]].x[1]
    y2 = ap.nd[ap.id[2]].x[2]
    return (y2-y1)/L,(x1-x2)/L,0.
end
function get_normal(ap1::AbstractSeg,ap2::AbstractPoi)
    id₁ = findfirst(x -> x == ap2.id[1], ap1.id)
    return (id₁ == 1 ? (-1.,0.,0.) : (1.,0.,0.))
end

# --------------- Seg2 ---------------
struct Seg2 <: AbstractSeg
    nd::Vector{PhysicalNode}
    id::Vector{Int}
    qw::Vector{ParametricNode}
    nm::Float64
end

# constructions of Seg2
function Seg2(nd::Vector{PhysicalNode},id::Vector{Int};qw::Symbol=:SegGI2)
    L = norm(nd[id[2]].x - nd[id[1]].x)
    qw = QuadratureRule[qw]
    return Seg2(nd,id,qw,L)
end
function Seg2(nd::Vector{PhysicalNode},ids::Vector{Vector{Int}};qw::Symbol=:SegGI2)
    return [Seg2(nd,id,qw=qw) for id in ids]
end

# actions of Seg2
get_shape_functions(::Seg2,ξ::AbstractVector{Float64},::Val{:∂1}) = SVector{2,Float64}((1.0-ξ[1])*0.5,(1.0+ξ[1])*0.5)
function get_shape_functions(ap::Seg2,::AbstractVector{Float64},::Val{:∂x})
    x1 = ap.nd[ap.id[1]].x[1]
    x2 = ap.nd[ap.id[2]].x[1]
    return SVector{2,Float64}(-1.0/(x2-x1),1.0/(x2-x1))
end
get_shape_functions(::Seg2,ξ::AbstractVector{Float64},::Val{:∂y}) = SVector{2,Float64}(0.,0.)
get_shape_functions(::Seg2,ξ::AbstractVector{Float64},::Val{:∂z}) = SVector{2,Float64}(0.,0.)

## AbstractTri
get_jacobe(ap::AbstractTri,ξ::AbstractVector{Float64}) = ap.nm
function get_coordinates(ap::AbstractTri,ξ::AbstractVector{Float64})
    return ξ[1]*ap.nd[ap.id[1]].x +
           ξ[2]*ap.nd[ap.id[2]].x +
           (1-ξ[1]-ξ[2])*ap.nd[ap.id[3]].x
end
function get_coordinates(ap1::AbstractTri,ap2::AbstractSeg,ξ::AbstractVector{Float64})
    id₁ = findfirst(x -> x == ap2.id[1], ap1.id)
    id₂ = findfirst(x -> x == ap2.id[2], ap1.id)
    return SVector{2,Float64}((1-ξ[1])/2*(id₁ == 1) + (1+ξ[1])/2*(id₂ == 1),
                              (1-ξ[1])/2*(id₁ == 2) + (1+ξ[1])/2*(id₂ == 2))
end
function get_normal(ap::AbstractTri)
    A = ap.nm
    x1 = ap.nd[ap.id[1]].x[1]
    y1 = ap.nd[ap.id[1]].x[2]
    z1 = ap.nd[ap.id[1]].x[3]
    x2 = ap.nd[ap.id[2]].x[1]
    y2 = ap.nd[ap.id[2]].x[2]
    z2 = ap.nd[ap.id[2]].x[3]
    x3 = ap.nd[ap.id[3]].x[1]
    y3 = ap.nd[ap.id[3]].x[2]
    z3 = ap.nd[ap.id[3]].x[3]
    Ax = 0.5*(y1*z2+y2*z3+y3*z1-y2*z1-y3x*z2-y1*z3)
    Ay = 0.5*(z1*x2+z2*x3+z3*x1-z2*x1-z3x*x2-z1*x3)
    Az = 0.5*(x1*y2+x2*y3+x3*y1-x2*y1-x3x*y2-x1*y3)
    return Ax/A,Ay/A,Az/A
end
function get_normal(ap1::AbstractTri,ap2::AbstractSeg)
    id₁ = findfirst(x -> x == ap2.id[1], ap1.id)
    id₂ = findfirst(x -> x == ap2.id[2], ap1.id)
    x1 = ap1.nd[ap1.id[id₁]].x[1]
    y1 = ap1.nd[ap1.id[id₁]].x[2]
    x2 = ap1.nd[ap1.id[id₂]].x[1]
    y2 = ap1.nd[ap1.id[id₂]].x[2]
    L = ap2.nm
    return (y2-y1)/L,(x1-x2)/L,0.
end

# --------------- Tri3 ---------------
# Constant strain triangular Approximator (CST)
struct Tri3 <: AbstractTri
    nd :: Vector{PhysicalNode}
    id :: Vector{Int}
    qw::Vector{ParametricNode}
    nm :: Float64
end

# constructions
function Tri3(x::Vector{PhysicalNode},id::Vector{Int};qw::Symbol=:TriGI3)
    x1 = x[id[1]].x[1]
    y1 = x[id[1]].x[2]
    z1 = x[id[1]].x[3]
    x2 = x[id[2]].x[1]
    y2 = x[id[2]].x[2]
    z2 = x[id[2]].x[3]
    x3 = x[id[3]].x[1]
    y3 = x[id[3]].x[2]
    z3 = x[id[3]].x[3]
    Ax = 0.5*(y1*z2+y2*z3+y3*z1-y2*z1-y3*z2-y1*z3)
    Ay = 0.5*(z1*x2+z2*x3+z3*x1-z2*x1-z3*x2-z1*x3)
    Az = 0.5*(x1*y2+x2*y3+x3*y1-x2*y1-x3*y2-x1*y3)
    A = (Ax^2 + Ay^2 + Az^2)^0.5
    qw = QuadratureRule[qw]
    return Tri3(x,id,qw,A)
end

# actions
get_shape_functions(ap::Tri3,ξ::AbstractVector{Float64},::Val{:∂1}) = SVector{3,Float64}(ξ[1],ξ[2],1-ξ[1]-ξ[2])
function get_shape_functions(ap::Tri3,ξ::AbstractVector{Float64},::Val{:∂x})
    y1 = ap.nd[ap.id[1]].x[2]
    y2 = ap.nd[ap.id[2]].x[2]
    y3 = ap.nd[ap.id[3]].x[2]
    A = ap.nm
    return SVector{3,Float64}((y2-y3)/(2A),(y3-y1)/(2A),(y1-y2)/(2A))
end
function get_shape_functions(ap::Tri3,ξ::AbstractVector{Float64},::Val{:∂y})
    x1 = ap.nd[ap.id[1]].x[1]
    x2 = ap.nd[ap.id[2]].x[1]
    x3 = ap.nd[ap.id[3]].x[1]
    A = ap.nm
    return SVector{3,Float64}((x3-x2)/(2A),(x1-x3)/(2A),(x2-x1)/(2A))
end
get_shape_functions(ap::Tri3,ξ::AbstractVector{Float64},::Val{:∂z}) = SVector{3,Float64}(0.,0.,0.)

## AbstractQuad
function get_jacobe(ap::AbstractQuad,ξ::AbstractVector{Float64})
    J₁₁,J₂₁,J₁₂,J₂₂ = get_jacobe_matrix(ap,ξ)
    return J₁₁*J₂₂-J₂₁*J₁₂
end
function get_jacobe_matrix(ap::AbstractQuad,ξ::AbstractVector{Float64})
    x₁ = ap.nd[ap.id[1]].x[1]
    x₂ = ap.nd[ap.id[2]].x[1]
    x₃ = ap.nd[ap.id[3]].x[1]
    x₄ = ap.nd[ap.id[4]].x[1]
    y₁ = ap.nd[ap.id[1]].x[2]
    y₂ = ap.nd[ap.id[2]].x[2]
    y₃ = ap.nd[ap.id[3]].x[2]
    y₄ = ap.nd[ap.id[4]].x[2]
    ∂N₁∂ξ,∂N₂∂ξ,∂N₃∂ξ,∂N₄∂ξ = get_shape_functions(ap,ξ,Val(:∂ξ))
    ∂N₁∂η,∂N₂∂η,∂N₃∂η,∂N₄∂η = get_shape_functions(ap,ξ,Val(:∂η))
    J₁₁ = ∂N₁∂ξ*x₁ + ∂N₂∂ξ*x₂ + ∂N₃∂ξ*x₃ + ∂N₄∂ξ*x₄
    J₁₂ = ∂N₁∂η*x₁ + ∂N₂∂η*x₂ + ∂N₃∂η*x₃ + ∂N₄∂η*x₄
    J₂₁ = ∂N₁∂ξ*y₁ + ∂N₂∂ξ*y₂ + ∂N₃∂ξ*y₃ + ∂N₄∂ξ*y₄
    J₂₂ = ∂N₁∂η*y₁ + ∂N₂∂η*y₂ + ∂N₃∂η*y₃ + ∂N₄∂η*y₄
    return J₁₁,J₂₁,J₁₂,J₂₂
end
function get_coordinates(ap::AbstractQuad,ξ::AbstractVector{Float64})
    N₁,N₂,N₃,N₄ = get_shape_functions(ap,ξ,Val(:∂1))
    return N₁*ap.nd[ap.id[1]] +
           N₂*ap.nd[ap.id[2]] +
           N₃*ap.nd[ap.id[3]] +
           N₄*ap.nd[ap.id[4]]
end
function get_coordinates(ap1::AbstractQuad,ap2::AbstractSeg,ξ::AbstractVector{Float64})
    id₁ = findfirst(x -> x == ap2.id[1], ap1.id)
    id₂ = findfirst(x -> x == ap2.id[2], ap1.id)
    return SVector{2,Float64}((1-ξ[1])/2*(id₁ == 1)*(id₂ == 2) + (1+ξ[1])/2*(id₁ == 3)*(id₂ == 4),
                              (1-ξ[1])/2*(id₁ == 2)*(id₂ == 3) + (1+ξ[1])/2*(id₁ == 4)*(id₂ == 1))
end

# --------------- Quad ---------------
mutable struct Quad <: AbstractQuad
    nd :: Vector{PhysicalNode}
    id :: Vector{Int}
    qw::Vector{ParametricNode}
end
# constructions
function Quad(x::Vector{PhysicalNode},id::Vector{Int};qw::Symbol=:QuadGI2)
    qw = QuadratureRule[qw]
    return Quad(x,id,qw)
end

# actions
function get_shape_functions(ap::Quad,ξ::AbstractVector{Float64},::Val{:∂1})
    N₁ = (1-ξ[1])*(1-ξ[2])/4
    N₂ = (1+ξ[1])*(1-ξ[2])/4
    N₃ = (1+ξ[1])*(1+ξ[2])/4
    N₄ = (1-ξ[1])*(1+ξ[2])/4
    return SVector{4,Float64}(N₁,N₂,N₃,N₄)
end
function get_shape_functions(ap::Quad,ξ::AbstractVector{Float64},::Val{:∂ξ})
    ∂N₁∂ξ = - (1-ξ[2])/4
    ∂N₂∂ξ =   (1-ξ[2])/4
    ∂N₃∂ξ =   (1+ξ[2])/4
    ∂N₄∂ξ = - (1+ξ[2])/4
    return SVector{4,Float64}(∂N₁∂ξ,∂N₂∂ξ,∂N₃∂ξ,∂N₄∂ξ)
end
function get_shape_functions(ap::Quad,ξ::AbstractVector{Float64},::Val{:∂η})
    ∂N₁∂η = - (1-ξ[1])/4
    ∂N₂∂η = - (1+ξ[1])/4
    ∂N₃∂η =   (1+ξ[1])/4
    ∂N₄∂η =   (1-ξ[1])/4
    return SVector{4,Float64}(∂N₁∂η,∂N₂∂η,∂N₃∂η,∂N₄∂η)
end
function get_shape_functions(ap::Quad,ξ::AbstractVector{Float64},::Val{:∂x},::Val{:∂y})
    x₁ = ap.nd[ap.id[1]].x[1]
    x₂ = ap.nd[ap.id[2]].x[1]
    x₃ = ap.nd[ap.id[3]].x[1]
    x₄ = ap.nd[ap.id[4]].x[1]
    y₁ = ap.nd[ap.id[1]].x[2]
    y₂ = ap.nd[ap.id[2]].x[2]
    y₃ = ap.nd[ap.id[3]].x[2]
    y₄ = ap.nd[ap.id[4]].x[2]
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
    return SVector{4,Float64}(∂N₁∂x,∂N₂∂x,∂N₃∂x,∂N₄∂x),SVector{4,Float64}(∂N₁∂y,∂N₂∂y,∂N₃∂y,∂N₄∂y)
end

get_shape_functions(ap::Quad,ξ::AbstractVector{Float64},::Val{:∂z}) = SVector{4,Float64}(0.,0.,0.,0.)
@inline get_shape_functions(ap::Quad,ξ::AbstractVector{Float64},::Val{:∂1},::Val{:∂x},::Val{:∂y}) = get_shape_functions(ap,ξ,Val(:∂1)),get_shape_functions(ap,ξ,Val(:∂x),Val(:∂y))...
@inline get_shape_functions(ap::Quad,ξ::AbstractVector{Float64},::Val{:∂1},::Val{:∂x},::Val{:∂y},::Val{:∂z}) = get_shape_functions(ap,ξ,Val(:∂1)),get_shape_functions(ap,ξ,Val(:∂x),Val(:∂y))...,get_shape_functions(ap,ξ,Val(:∂z))
