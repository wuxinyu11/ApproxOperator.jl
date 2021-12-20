## Actions for general functions
@inline get_global_indice(ap::Approximator,i::Int) = ap.id[i]
@inline get_number_of_indices(ap::Approximator) = length(ap.id)
@inline get_integration_points_and_weights(ap::Approximator) = ap.qw
@inline get_local_node(ap::Approximator,i::Int) = ap.nodes[ap.id[i]]
# @inline get_shape_functions(ap::Approximator,ξ::Union{Float64,AbstractVector},gs::Val...) = (get_shape_functions(ap,ξ,g) for g in gs)
@inline get_shape_functions(ap::Approximator,ξ::Union{Float64,AbstractVector},::Val{:∂1},::Val{:∂x}) = (get_shape_functions(ap,ξ,Val(:∂1)),get_shape_functions(ap,ξ,Val(:∂x)))
@inline get_shape_functions(ap::Approximator,ξ::Union{Float64,AbstractVector},::Val{:∂1},::Val{:∂x},::Val{:∂y}) = (get_shape_functions(ap,ξ,Val(:∂1)),get_shape_functions(ap,ξ,Val(:∂x)),get_shape_functions(ap,ξ,Val(:∂y)))
@inline get_shape_functions(ap::Approximator,ξ::Union{Float64,AbstractVector},::Val{:∂1},::Val{:∂x},::Val{:∂y},::Val{:∂z}) = (get_shape_functions(ap,ξ,Val(:∂1)),get_shape_functions(ap,ξ,Val(:∂x)),get_shape_functions(ap,ξ,Val(:∂y)),get_shape_functions(ap,ξ,Val(:∂z)))

function set_integration_rule!(ap::Approximator,qw::Symbol)
    ap.qw = QuadratureRule[qw]
end
function set_integration_rule!(aps::Vector{Approximator},qw::Symbol)
    for ap in aps
        set_integration_rule!(ap,qw)
    end
end
## AbstractPoi
get_jacobe(::AbstractPoi,::Float64) = 1.
get_coordinates(ap::AbstractPoi,::Float64) = 1.0*ap.nodes[ap.id[1]]

# --------------- Poi1 ---------------
struct Poi1 <: AbstractPoi
    nodes::Vector{Node}
    id::Int
    qw::Vector{Pair{Float64,Float64}}
end

# constructions of Poi1
Poi1(nodes::Vector{Node},id::Int;qw::Symbol=:PoiGI1) = Poi1(nodes,id,QuadratureRule[qw])
function Poi1(nodes::Vector{Node},ids::Vector{Int};qw::Symbol=:PoiGI1)
    return [Poi1(nodes,id,qw=qw) for id in ids]
end

# actions of Poi1
get_shape_functions(::Poi1,::Float64,::Val{:∂1}) = 1.
get_shape_functions(::Poi1,::Float64,::Val{:∂x}) = 1.


## AbstractSeg
get_number_of_dimensions(::AbstractSeg) = 1
get_jacobe(ap::AbstractSeg,ξ::Float64) = ap.norm/2
function get_coordinates(ap::AbstractSeg,ξ::Float64)
    N1 = (1.0-ξ)*0.5
    N2 = (1.0+ξ)*0.5
    return N1*ap.nodes[ap.id[1]] + N2*ap.nodes[ap.id[2]]
end
function get_coordinates(ap1::AbstractSeg,ap2::AbstractPoi,::Float64)
    id₁ = findfirst(x -> x == ap2.id[1], ap1.id)
    return (id₁ == 1 ? -1. : 1.)
end
function get_normal(ap::AbstractSeg)
    L = ap.norm
    x1 = ap.nodes[ap.id[1]].x
    y1 = ap.nodes[ap.id[1]].y
    x2 = ap.nodes[ap.id[2]].x
    y2 = ap.nodes[ap.id[2]].y
    return (y2-y1)/L,(x1-x2)/L,0.
end
function get_normal(ap1::AbstractSeg,ap2::AbstractPoi)
    id₁ = findfirst(x -> x == ap2.id[1], ap1.id)
    return (id₁ == 1 ? (-1.,0.,0.) : (1.,0.,0.))
end

# --------------- Seg2 ---------------
struct Seg2 <: AbstractSeg
    nodes::Vector{Node}
    id::Vector{Int}
    qw::Vector{Pair{Float64,Float64}}
    norm :: Float64
end

# constructions of Seg2
function Seg2(nodes::Vector{Node},id::Vector{Int};qw::Symbol=:SegGI2)
    L = norm(nodes[id[2]] - nodes[id[1]])
    qw = QuadratureRule[qw]
    return Seg2(nodes,id,qw,L)
end
function Seg2(nodes::Vector{Node},ids::Vector{Vector{Int}};qw::Symbol=:SegGI2)
    return [Seg2(nodes,id,qw=qw) for id in ids]
end

# actions of Seg2
get_shape_functions(::Seg2,ξ::Float64,::Val{:∂1}) = SVector{2,Float64}((1.0-ξ)*0.5,(1.0+ξ)*0.5)
function get_shape_functions(ap::Seg2,ξ::Float64,::Val{:∂x})
    x1 = ap.nodes[ap.id[1]].x
    x2 = ap.nodes[ap.id[2]].x
    return SVector{2,Float64}(-1.0/(x2-x1),1.0/(x2-x1))
end
get_shape_functions(::Seg2,ξ::Float64,::Val{:∂y}) = SVector{2,Float64}(0.,0.)
get_shape_functions(::Seg2,ξ::Float64,::Val{:∂z}) = SVector{2,Float64}(0.,0.)

## AbstractTri
get_number_of_dimensions(::AbstractTri) = 2
get_jacobe(ap::AbstractTri,ξ::Vector{Float64}) = ap.norm
function get_coordinates(ap::AbstractTri,ξ::AbstractVector{Float64})
    return ξ[1]*ap.nodes[ap.id[1]] +
           ξ[2]*ap.nodes[ap.id[2]] +
           (1-ξ[1]-ξ[2])*ap.nodes[ap.id[3]]
end
function get_coordinates(ap1::AbstractTri,ap2::AbstractSeg,ξ::Float64)
    id₁ = findfirst(x -> x == ap2.id[1], ap1.id)
    id₂ = findfirst(x -> x == ap2.id[2], ap1.id)
    return SVector{2,Float64}((1-ξ)/2*(id₁ == 1) + (1+ξ)/2*(id₂ == 1),
                              (1-ξ)/2*(id₁ == 2) + (1+ξ)/2*(id₂ == 2))
end
function get_normal(ap::AbstractTri)
    A = ap.norm
    x1 = ap.nodes[ap.id[1]].x
    y1 = ap.nodes[ap.id[1]].y
    z1 = ap.nodes[ap.id[1]].z
    x2 = ap.nodes[ap.id[2]].x
    y2 = ap.nodes[ap.id[2]].y
    z2 = ap.nodes[ap.id[2]].z
    x3 = ap.nodes[ap.id[3]].x
    y3 = ap.nodes[ap.id[3]].y
    z3 = ap.nodes[ap.id[3]].z
    Ax = 0.5*(y1*z2+y2*z3+y3*z1-y2*z1-y3x*z2-y1*z3)
    Ay = 0.5*(z1*x2+z2*x3+z3*x1-z2*x1-z3x*x2-z1*x3)
    Az = 0.5*(x1*y2+x2*y3+x3*y1-x2*y1-x3x*y2-x1*y3)
    return Ax/A,Ay/A,Az/A
end
function get_normal(ap1::AbstractTri,ap2::AbstractSeg)
    id₁ = findfirst(x -> x == ap2.id[1], ap1.id)
    id₂ = findfirst(x -> x == ap2.id[2], ap1.id)
    x1 = ap1.nodes[ap1.id[id₁]].x
    y1 = ap1.nodes[ap1.id[id₁]].y
    x2 = ap1.nodes[ap1.id[id₂]].x
    y2 = ap1.nodes[ap1.id[id₂]].y
    L = ap2.norm
    return (y2-y1)/L,(x1-x2)/L,0.
end

# --------------- Tri3 ---------------
# Constant strain triangular Approximator (CST)
struct Tri3 <: AbstractTri
    nodes :: Vector{Node}
    id :: Vector{Int}
    qw::Vector{Pair{Vector{Float64},Float64}}
    norm :: Float64
end

# constructions
function Tri3(x::Vector{Node},id::Vector{Int};qw::Symbol=:TriGI3)
    x1 = x[id[1]].x
    y1 = x[id[1]].y
    z1 = x[id[1]].z
    x2 = x[id[2]].x
    y2 = x[id[2]].y
    z2 = x[id[2]].z
    x3 = x[id[3]].x
    y3 = x[id[3]].y
    z3 = x[id[3]].z
    Ax = 0.5*(y1*z2+y2*z3+y3*z1-y2*z1-y3*z2-y1*z3)
    Ay = 0.5*(z1*x2+z2*x3+z3*x1-z2*x1-z3*x2-z1*x3)
    Az = 0.5*(x1*y2+x2*y3+x3*y1-x2*y1-x3*y2-x1*y3)
    A = (Ax^2 + Ay^2 + Az^2)^0.5
    qw = QuadratureRule[qw]
    return Tri3(x,id,qw,A)
end
function Tri3(x::Vector{Node},ids::Vector{Vector{Int}};qw::Symbol=:TriGI3)
    return [Tri3(x,id,qw=qw) for id in ids]
end

# actions
get_shape_functions(ap::Tri3,ξ::AbstractVector{Float64},::Val{:∂1}) = SVector{3,Float64}(ξ[1],ξ[2],1-ξ[1]-ξ[2])
function get_shape_functions(ap::Tri3,ξ::AbstractVector{Float64},::Val{:∂x})
    y1 = ap.nodes[ap.id[1]].y
    y2 = ap.nodes[ap.id[2]].y
    y3 = ap.nodes[ap.id[3]].y
    A = ap.norm
    return SVector{3,Float64}((y2-y3)/(2A),(y3-y1)/(2A),(y1-y2)/(2A))
end
function get_shape_functions(ap::Tri3,ξ::AbstractVector{Float64},::Val{:∂y})
    x1 = ap.nodes[ap.id[1]].x
    x2 = ap.nodes[ap.id[2]].x
    x3 = ap.nodes[ap.id[3]].x
    A = ap.norm
    return SVector{3,Float64}((x3-x2)/(2A),(x1-x3)/(2A),(x2-x1)/(2A))
end
get_shape_functions(ap::Tri3,ξ::AbstractVector{Float64},::Val{:∂z}) = SVector{3,Float64}(0.,0.,0.)

## AbstractQuad
get_number_of_dimensions(::AbstractQuad) = 2
function get_jacobe(ap::AbstractQuad,ξ::Vector{Float64})
    J₁₁,J₂₁,J₁₂,J₂₂ = get_jacobe_matrix(ap,ξ)
    return J₁₁*J₂₂-J₂₁*J₁₂
end
function get_jacobe_matrix(ap::AbstractQuad,ξ::Vector{Float64})
    x₁ = ap.nodes[ap.id[1]].x
    x₂ = ap.nodes[ap.id[2]].x
    x₃ = ap.nodes[ap.id[3]].x
    x₄ = ap.nodes[ap.id[4]].x
    y₁ = ap.nodes[ap.id[1]].y
    y₂ = ap.nodes[ap.id[2]].y
    y₃ = ap.nodes[ap.id[3]].y
    y₄ = ap.nodes[ap.id[4]].y
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
    return N₁*ap.nodes[ap.id[1]] +
           N₂*ap.nodes[ap.id[2]] +
           N₃*ap.nodes[ap.id[3]] +
           N₄*ap.nodes[ap.id[4]]
end
function get_coordinates(ap1::AbstractQuad,ap2::AbstractSeg,ξ::Float64)
    id₁ = findfirst(x -> x == ap2.id[1], ap1.id)
    id₂ = findfirst(x -> x == ap2.id[2], ap1.id)
    return SVector{2,Float64}((1-ξ)/2*(id₁ == 1)*(id₂ == 2) + (1+ξ)/2*(id₁ == 3)*(id₂ == 4),
                              (1-ξ)/2*(id₁ == 2)*(id₂ == 3) + (1+ξ)/2*(id₁ == 4)*(id₂ == 1))
end

# --------------- Quad ---------------
mutable struct Quad <: AbstractQuad
    nodes :: Vector{Node}
    id :: Vector{Int}
    qw::Vector{Pair{Vector{Float64},Float64}}
end
# constructions
function Quad(x::Vector{Node},id::Vector{Int};qw::Symbol=:QuadGI2)
    qw = QuadratureRule[qw]
    return Quad(x,id,qw)
end

function Quad(x::Vector{Node},ids::Vector{Vector{Int}};qw::Symbol=:QuadGI2)
    return [Quad(x,id,qw=qw) for id in ids]
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
    x₁ = ap.nodes[ap.id[1]].x
    x₂ = ap.nodes[ap.id[2]].x
    x₃ = ap.nodes[ap.id[3]].x
    x₄ = ap.nodes[ap.id[4]].x
    y₁ = ap.nodes[ap.id[1]].y
    y₂ = ap.nodes[ap.id[2]].y
    y₃ = ap.nodes[ap.id[3]].y
    y₄ = ap.nodes[ap.id[4]].y
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

get_shape_functions(ap::Quad,ξ::AbstractVector{Float64},::Val{:∂1},::Val{:∂x},::Val{:∂y}) = get_shape_functions(ap,ξ,Val(:∂1)),get_shape_functions(ap,ξ,Val(:∂x),Val(:∂y))...
get_shape_functions(ap::Quad,ξ::AbstractVector{Float64},::Val{:∂1},::Val{:∂x},::Val{:∂y},::Val{:∂z}) = get_shape_functions(ap,ξ,Val(:∂1)),get_shape_functions(ap,ξ,Val(:∂x),Val(:∂y))...,get_shape_functions(ap,ξ,Val(:∂z))

## Meshfree

## BasisFunction
# ------------ Linear1D ---------------
struct Linear1D <: BasisFunction
    𝗠::Dict{Symbol,SymMat}
end
# constructions of BasisFunction
function Linear1D(gs::Symbol...)
    𝗠 = Dict{Symbol,SymMat}()
    for g in gs
        push!(𝗠,g=>SymMat(2))
    end
    return Linear1D(𝗠)
end

# actions of BasisFunction
@inline get_basis_function(::Linear1D,x::AbstractVector,::Val{:∂1}) = SVector{2,Float64}(1.,x[1])
@inline get_basis_function(::Linear1D,::AbstractVector,::Val{:∂x}) = SVector{2,Float64}(0.,1.)
@inline get_basis_function(::Linear1D,::AbstractVector,::Val{:∂y}) = SVector{2,Float64}(0.,0.)
@inline get_basis_function(::Linear1D,::AbstractVector,::Val{:∂z}) = SVector{2,Float64}(0.,0.)

# ------------ Quadaratic1D ---------------
struct Quadratic1D <: BasisFunction
    𝗠::Dict{Symbol,SymMat}
end

# constructions of BasisFunction
function Quadratic1D(gs::Symbol...)
    𝗠 = Dict{Symbol,SymMat}()
    for g in gs
        push!(𝗠,g=>SymMat(3))
    end
    return Quadratic1D(𝗠)
end

# actions of BasisFunction
@inline get_basis_function(::Quadratic1D,x::AbstractVector,::Val{:∂1}) = SVector{3,Float64}(1.,x[1],x[1]^2)
@inline get_basis_function(::Quadratic1D,x::AbstractVector,::Val{:∂x}) = SVector{3,Float64}(0.,1.,2*x[1])
@inline get_basis_function(::Quadratic1D,x::AbstractVector,::Val{:∂y}) = SVector{3,Float64}(0.,0.,0.)
@inline get_basis_function(::Quadratic1D,x::AbstractVector,::Val{:∂z}) = SVector{3,Float64}(0.,0.,0.)
@inline get_basis_function(::Quadratic1D,x::AbstractVector,::Val{:∂x²}) = SVector{3,Float64}(0.,0.,2.)

# ------------ Cubic1D ---------------
struct Cubic1D <: BasisFunction
    𝗠::Dict{Symbol,SymMat}
end

# constructions of BasisFunction
function Cubic1D(gs::Symbol...)
    𝗠 = Dict{Symbol,SymMat}()
    for g in gs
        push!(𝗠,g=>SymMat(4))
    end
    return Cubic1D(𝗠)
end

# actions of BasisFunction
@inline get_basis_function(::Cubic1D,x::AbstractVector,::Val{:∂1}) = SVector{4,Float64}(1.,x[1],x[1]^2,x[1]^3)
@inline get_basis_function(::Cubic1D,x::AbstractVector,::Val{:∂x}) = SVector{4,Float64}(0.,1.,2*x[1],3*x[1]^2)
@inline get_basis_function(::Cubic1D,x::AbstractVector,::Val{:∂y}) = SVector{4,Float64}(0.,0.,0.,0.)
@inline get_basis_function(::Cubic1D,x::AbstractVector,::Val{:∂z}) = SVector{4,Float64}(0.,0.,0.,0.)
@inline get_basis_function(::Cubic1D,x::AbstractVector,::Val{:∂x²}) = SVector{4,Float64}(0.,0.,2.,6*x[1])

# ------------ Linear2D ---------------
struct Linear2D <: BasisFunction
    𝗠::Dict{Symbol,SymMat}
end
# constructions of BasisFunction
function Linear2D(gs::Symbol...)
    𝗠 = Dict{Symbol,SymMat}()
    for g in gs
        push!(𝗠,g=>SymMat(3))
    end
    return Linear2D(𝗠)
end

# actions of BasisFunction
@inline get_basis_function(::Linear2D,x::AbstractVector,::Val{:∂1}) = SVector{3,Float64}(1.,x[1],x[2])
@inline get_basis_function(::Linear2D,::AbstractVector,::Val{:∂x}) = SVector{3,Float64}(0.,1.,0.)
@inline get_basis_function(::Linear2D,::AbstractVector,::Val{:∂y}) = SVector{3,Float64}(0.,0.,1.)
@inline get_basis_function(::Linear2D,::AbstractVector,::Val{:∂z}) = SVector{3,Float64}(0.,0.,0.)

# ------------ Quadratic2D ---------------
struct Quadratic2D <: BasisFunction
    𝗠::Dict{Symbol,SymMat}
end
# constructions of BasisFunction
function Quadratic2D(gs::Symbol...)
    𝗠 = Dict{Symbol,SymMat}()
    for g in gs
        push!(𝗠,g=>SymMat(6))
    end
    return Quadratic2D(𝗠)
end

# actions of BasisFunction
@inline get_basis_function(::Quadratic2D,x::AbstractVector,::Val{:∂1}) = SVector{6,Float64}(1.,x[1],x[2],x[1]^2,x[1]*x[2],x[2]^2)
@inline get_basis_function(::Quadratic2D,x::AbstractVector,::Val{:∂x}) = SVector{6,Float64}(0.,1.,0.,2*x[1],x[2],0.)
@inline get_basis_function(::Quadratic2D,x::AbstractVector,::Val{:∂y}) = SVector{6,Float64}(0.,0.,1.,0.,x[1],2*x[2])
@inline get_basis_function(::Quadratic2D,::AbstractVector,::Val{:∂z}) = SVector{6,Float64}(0.,0.,0.,0.,0.,0.)

# ------------ Cubic2D ---------------
struct Cubic2D <: BasisFunction
    𝗠::Dict{Symbol,SymMat}
end
# constructions of BasisFunction
function Cubic2D(gs::Symbol...)
    𝗠 = Dict{Symbol,SymMat}()
    for g in gs
        push!(𝗠,g=>SymMat(6))
    end
    return Cubic2D(𝗠)
end

# actions of BasisFunction
@inline get_basis_function(::Cubic2D,x::AbstractVector,::Val{:∂1}) =
SVector{10,Float64}(
    1., x[1], x[2], x[1]^2, x[1]*x[2], x[2]^2, x[1]^3, x[1]^2*x[2], x[1]*x[2]^2, x[2]^3
)
@inline get_basis_function(::Cubic2D,x::AbstractVector,::Val{:∂x}) =
SVector{10,Float64}(
    0., 1., 0., 2*x[1], x[2], 0., 3*x[1]^2, 2*x[1]*x[2], x[2]^2, 0.
)
@inline get_basis_function(::Cubic2D,x::AbstractVector,::Val{:∂y}) =
SVector{10,Float64}(
    0., 0., 1., 0., x[1], 2*x[2], 0., x[1]^2, 2*x[1]*x[2], 3*x[2]^2
)
@inline get_basis_function(::Cubic2D,::AbstractVector,::Val{:∂z}) =
SVector{10,Float64}(
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
)

# --------------- TensorProductKernel ---------------
struct TensorProductKernel <: KernelFunction
    support_size::Vector{Float64}
    kernel_type::Symbol
    𝝭::Dict{Symbol,Vector{Float64}}
end

# constructions of TensorProductKernel
function TensorProductKernel(gs::Symbol...;ss::Vector{Float64}=[1.,1.,1.],nm::Int=10,kt::Symbol=:CubicSpline)
    𝝭 = Dict{Symbol,Vector{Float64}}()
    for g in gs
        push!(𝝭,g=>zeros(nm))
    end
    return TensorProductKernel(ss,kt,𝝭)
end

# actions of TensorProductKernel
function get_kernel_function(kf::TensorProductKernel,Δx::AbstractVector,::Val{:∂1})
    sᵢ = kf.support_size
    kt = Val(kf.kernel_type)
    rx = abs(Δx[1])/sᵢ[1]
    ry = abs(Δx[2])/sᵢ[2]
    rz = abs(Δx[3])/sᵢ[3]
    wx = get_kernel(kt,rx,Val(:∂1))
    wy = get_kernel(kt,ry,Val(:∂1))
    wz = get_kernel(kt,rz,Val(:∂1))
    return wx*wy*wz
end

function get_kernel_function(kf::TensorProductKernel,Δx::AbstractVector,::Val{:∂1},::Val{:∂x})
    sᵢ = kf.support_size
    kt = Val(kf.kernel_type)
    rx = abs(Δx[1])/sᵢ[1]
    ∂rx = sign(Δx[1])/sᵢ[1]
    wx = get_kernel(kt,rx,Val(:∂1))
    ∂wx = get_kernel(kt,rx,Val(:∂r))*∂rx
    return wx, ∂wx
end

function get_kernel_function(kf::TensorProductKernel,Δx::AbstractVector,::Val{:∂1},::Val{:∂x},::Val{:∂y})
    sᵢ = kf.support_size
    kt = Val(kf.kernel_type)
    rx = abs(Δx[1])/sᵢ[1]
    ry = abs(Δx[2])/sᵢ[2]
    ∂rx = sign(Δx[1])/sᵢ[1]
    ∂ry = sign(Δx[2])/sᵢ[2]
    wx = get_kernel(kt,rx,Val(:∂1))
    wy = get_kernel(kt,ry,Val(:∂1))
    ∂wx = get_kernel(kt,rx,Val(:∂r))*∂rx
    ∂wy = get_kernel(kt,ry,Val(:∂r))*∂ry
    return wx*wy, ∂wx*wy, wx*∂wy
end

function get_kernel_function(kf::TensorProductKernel,Δx::AbstractVector,::Val{:∂1},::Val{:∂x},::Val{:∂y},::Val{:∂z})
    sᵢ = kf.support_size
    kt = Val(kf.kernel_type)
    rx = abs(Δx[1])/sᵢ[1]
    ry = abs(Δx[2])/sᵢ[2]
    rz = abs(Δx[3])/sᵢ[3]
    ∂rx = sign(Δx[1])/sᵢ[1]
    ∂ry = sign(Δx[2])/sᵢ[2]
    ∂rz = sign(Δx[3])/sᵢ[3]
    wx = get_kernel(kt,rx,Val(:∂1))
    wy = get_kernel(kt,ry,Val(:∂1))
    wz = get_kernel(kt,rz,Val(:∂1))
    ∂wx = get_kernel(kt,rx,Val(:∂r))*∂rx
    ∂wy = get_kernel(kt,ry,Val(:∂r))*∂ry
    ∂wz = get_kernel(kt,rz,Val(:∂r))*∂rz
    return wx*wy*wz, ∂wx*wy*wz, wx*∂wy*wz, wx*wy*∂wz
end

# function get_kernel_function(kf::TensorProductKernel,Δx::SVector{3,Float64})
#     sᵢ = kf.support_size
#     kt = kf.kernel_type
#     rx = abs(Δx[1])/sᵢ[1]
#     ry = abs(Δx[2])/sᵢ[2]
#     rz = abs(Δx[3])/sᵢ[3]
#     ∂rx = sign(Δx[1])/sᵢ[1]
#     ∂ry = sign(Δx[2])/sᵢ[2]
#     ∂rz = sign(Δx[3])/sᵢ[3]
#     wx = get_kernel(kt,rx)
#     wy = get_kernel(kt,ry)
#     wz = get_kernel(kt,rz)
#     ∂wx = get_gradient_of_kernel(kt,rx)*∂rx
#     ∂wy = get_gradient_of_kernel(kt,ry)*∂ry
#     ∂wz = get_gradient_of_kernel(kt,rz)*∂rz
#     ∂²wx = get_2nd_gradient_of_kernel(kt,rx)*∂rx^2
#     ∂²wy = get_2nd_gradient_of_kernel(kt,ry)*∂ry^2
#     ∂²wz = get_2nd_gradient_of_kernel(kt,rz)*∂rz^2
#     return SVector{6,Float64}(∂²wx*wy*wz,
#                               ∂wx*∂wy*wz,
#                               wx*∂²wy*wz,
#                               ∂wx*wy*∂wz,
#                               wx*∂wy*∂wz,
#                               wx*wy*∂²wz)
# end
# ----------------- CircularKernel ---------------
struct CircularKernel <: KernelFunction
    support_size::Float64
    kernel_type::Symbol
    𝝭::Dict{Symbol,Vector{Float64}}
end

# --------------- Kernel ---------------
function get_kernel(::Val{:CubicSpline},r::Float64,::Val{:∂1})
    if r > 1.0
        return 0.0
    elseif r <= 0.5
        return 2/3 - 4*r^2 +  4*r^3
    else
        return 4/3 - 4*r + 4*r^2 - 4*r^3/3
    end
end

function get_kernel(::Val{:CubicSpline},r::Float64,::Val{:∂r})
    if r > 1.0
        return 0.0
    elseif r <= 0.5
        return - 8*r + 12*r^2
    else
        return - 4   + 8*r - 4*r^2
    end
end

function get_kernel(::Val{:CubicSpline},r::Float64,::Val{:∂r²})
    if r > 1.0
        return 0.0
    elseif r <= 0.5
        return - 8 + 24*r
    else
        return   8 - 8*r
    end
end

## shape function
struct RKShape<:ShapeFunction
    𝝭::Dict{Symbol,Dict{AbstractVector,SparseVector}}
end
## Meshfree
# -------------- PoiM ---------------
struct PoiM{B<:BasisFunction,K<:KernelFunction} <: AbstractPoi
    nodes::Vector{Node}
    id::Vector{Int}
    qw::Vector{Pair{Float64,Float64}}
    bf::B
    kf::K
end

# constructions of PoiM
function PoiM(nodes::Vector{Node},id::Vector{Int};bf::BasisFunction=Linear1D(),kf::KernelFunction=TensorProductKernel())
    return PoiM(nodes,id,bf=bf,kf=kf)
end
function PoiM(nodes::Vector{Node},id::Int;bf::BasisFunction=Linear1D(),kf::KernelFunction=TensorProductKernel(),sp::Union{SpatialPartition,Nothing}=nothing)
    if sp ≠ nothing
        id = union!([id],collect(sp(nodes[id])))
    end
    qw = QuadratureRule[:PoiGI1]
    return PoiM(nodes,id,qw,bf,kf)
end

# -------------- SegM ---------------
struct SegM{B<:BasisFunction,K<:KernelFunction} <: AbstractSeg
    nodes :: Vector{Node}
    id :: Vector{Int}
    qw::Vector{Pair{Float64,Float64}}
    norm::Float64
    bf::B
    kf::K
end
function SegM(nodes::Vector{Node},ids::Vector{Vector{Int}};qw::Symbol=:SegGI2,bf::BasisFunction=Linear1D(),kf::KernelFunction=TensorProductKernel(),sp::Union{SpatialPartition,Nothing}=nothing)
    return [SegM(nodes,id,qw=qw,bf=bf,kf=kf,sp=sp) for id in ids]
end
function SegM(nodes::Vector{Node},id::Vector{Int};qw::Symbol=:SegGI2,bf::BasisFunction=Linear1D(),kf::KernelFunction=TensorProductKernel(),sp::Union{SpatialPartition,Nothing}=nothing)
    if sp ≠ nothing
        id = union!(id,collect(sp(nodes[id])))
    end
    L = norm(nodes[id[2]] - nodes[id[1]])
    qw = QuadratureRule[qw]
    return SegM(nodes,id,qw,L,bf,kf)
end

# --------------- TriM ---------------
struct TriM{B<:BasisFunction,K<:KernelFunction} <: AbstractTri
    nodes :: Vector{Node}
    id :: Vector{Int}
    qw::Vector{Pair{Vector{Float64},Float64}}
    norm :: Float64
    bf:: B
    kf:: K
end

# constructions
function TriM(x::Vector{Node},ids::Vector{Vector{Int}};qw::Symbol=:TriGI3,bf::BasisFunction=Linear2D(),kf::KernelFunction=TensorProductKernel(),sp::Union{SpatialPartition,Nothing}=nothing)
    return [TriM(x,id,qw=qw,bf=bf,kf=kf,sp=sp) for id in ids]
end
function TriM(x::Vector{Node},id::Vector{Int};qw::Symbol=:TriGI3,bf::BasisFunction=Linear2D(),kf::KernelFunction=TensorProductKernel(),sp::Union{SpatialPartition,Nothing}=nothing)
    if sp ≠ nothing
        id = union!(id,collect(sp(x[id])))
    end
    x1 = x[id[1]].x
    y1 = x[id[1]].y
    z1 = x[id[1]].z
    x2 = x[id[2]].x
    y2 = x[id[2]].y
    z2 = x[id[2]].z
    x3 = x[id[3]].x
    y3 = x[id[3]].y
    z3 = x[id[3]].z
    Ax = 0.5*(y1*z2+y2*z3+y3*z1-y2*z1-y3*z2-y1*z3)
    Ay = 0.5*(z1*x2+z2*x3+z3*x1-z2*x1-z3*x2-z1*x3)
    Az = 0.5*(x1*y2+x2*y3+x3*y1-x2*y1-x3*y2-x1*y3)
    A = (Ax^2 + Ay^2 + Az^2)^0.5
    qw = QuadratureRule[qw]
    return TriM(x,id,qw,A,bf,kf)
end
# -------------- ReproducingKernel ---------------
# actions of ReproducingKernel
ReproducingKernel = Union{SegM{B,K},PoiM{B,K},TriM{B,K}} where {B,K}
function get_shape_functions(ap::ReproducingKernel,ξ::Union{Float64,AbstractVector{Float64}},::Val{:∂1})
    x = get_coordinates(ap,ξ)
    p₀ᵀ𝗠⁻¹ = cal_moment_matrix!(ap,x,Val(:∂1))
    𝝭 = get_shape_function(ap,:∂1)
    for i in 1:get_number_of_indices(ap)
        xᵢ = get_local_node(ap,i)
        Δx = x - xᵢ
        p = get_basis_function(ap.bf,Δx,Val(:∂1))
        w = get_kernel_function(ap.kf,Δx,Val(:∂1))
        𝝭[i] = p₀ᵀ𝗠⁻¹*p*w
    end
    return 𝝭
end

function get_shape_functions(ap::ReproducingKernel,ξ::Union{Float64,AbstractVector{Float64}},::Val{:∂1},::Val{:∂x})
    x = get_coordinates(ap,ξ)
    p₀ᵀ𝗠⁻¹, p₀ᵀ∂𝗠⁻¹∂x = cal_moment_matrix!(ap,x,Val(:∂1),Val(:∂x))
    # 𝝭, ∂𝝭∂x = get_shape_function(ap,:∂1,:∂x)
    𝝭 = get_shape_function(ap,:∂1)
    ∂𝝭∂x = get_shape_function(ap,:∂x)
    for i in 1:get_number_of_indices(ap)
        xᵢ = get_local_node(ap,i)
        Δx = x - xᵢ
        p = get_basis_function(ap.bf,Δx,Val(:∂1))
        ∂p∂x = get_basis_function(ap.bf,Δx,Val(:∂x))
        w, ∂w∂x = get_kernel_function(ap.kf,Δx,Val(:∂1),Val(:∂x))
        𝝭[i] = p₀ᵀ𝗠⁻¹*p*w
        ∂𝝭∂x[i] = p₀ᵀ∂𝗠⁻¹∂x*p*w + p₀ᵀ𝗠⁻¹*∂p∂x*w + p₀ᵀ𝗠⁻¹*p*∂w∂x
    end
    return 𝝭, ∂𝝭∂x
end

function get_shape_functions(ap::ReproducingKernel,ξ::AbstractVector{Float64},::Val{:∂1},::Val{:∂x},::Val{:∂y})
    x = get_coordinates(ap,ξ)
    p₀ᵀ𝗠⁻¹, p₀ᵀ∂𝗠⁻¹∂x, p₀ᵀ∂𝗠⁻¹∂y = cal_moment_matrix!(ap,x,Val(:∂1),Val(:∂x),Val(:∂y))
    # 𝝭, ∂𝝭∂x, ∂𝝭∂y = get_shape_function(ap,:∂1,:∂x,:∂y)
    𝝭 = get_shape_function(ap,:∂1)
    ∂𝝭∂x = get_shape_function(ap,:∂x)
    ∂𝝭∂y = get_shape_function(ap,:∂y)
    for i in 1:get_number_of_indices(ap)
        xᵢ = get_local_node(ap,i)
        Δx = x - xᵢ
        # p, ∂p∂x, ∂p∂y = get_basis_function(ap,Δx,Val(:∂1),Val(:∂x),Val(:∂y))
        p = get_basis_function(ap.bf,Δx,Val(:∂1))
        ∂p∂x = get_basis_function(ap.bf,Δx,Val(:∂x))
        ∂p∂y = get_basis_function(ap.bf,Δx,Val(:∂y))
        w, ∂w∂x, ∂w∂y = get_kernel_function(ap.kf,Δx,Val(:∂1),Val(:∂x),Val(:∂y))
        𝝭[i] = p₀ᵀ𝗠⁻¹*p*w
        ∂𝝭∂x[i] = p₀ᵀ∂𝗠⁻¹∂x*p*w + p₀ᵀ𝗠⁻¹*∂p∂x*w + p₀ᵀ𝗠⁻¹*p*∂w∂x
        ∂𝝭∂y[i] = p₀ᵀ∂𝗠⁻¹∂y*p*w + p₀ᵀ𝗠⁻¹*∂p∂y*w + p₀ᵀ𝗠⁻¹*p*∂w∂y
    end
    return 𝝭, ∂𝝭∂x, ∂𝝭∂y
end

function get_shape_functions(ap::ReproducingKernel,ξ::AbstractVector{Float64},::Val{:∂1},::Val{:∂x},::Val{:∂y},::Val{:∂z})
    x = get_coordinates(ap,ξ)
    p₀ᵀ𝗠⁻¹, p₀ᵀ∂𝗠⁻¹∂x, p₀ᵀ∂𝗠⁻¹∂y, p₀ᵀ∂𝗠⁻¹∂z = cal_moment_matrix!(ap,x,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
    # 𝝭, ∂𝝭∂x, ∂𝝭∂y, ∂𝝭∂z = get_shape_function(ap,:∂1,:∂x,:∂y,:∂z)
    𝝭 = get_shape_function(ap,:∂1)
    ∂𝝭∂x = get_shape_function(ap,:∂x)
    ∂𝝭∂y = get_shape_function(ap,:∂y)
    ∂𝝭∂z = get_shape_function(ap,:∂z)
    for i in 1:get_number_of_indices(ap)
        xᵢ = get_local_node(ap,i)
        Δx = x - xᵢ
        # p, ∂p∂x, ∂p∂y, ∂p∂z = get_basis_function(ap,Δx,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
        p = get_basis_function(ap.bf,Δx,Val(:∂1))
        ∂p∂x = get_basis_function(ap.bf,Δx,Val(:∂x))
        ∂p∂y = get_basis_function(ap.bf,Δx,Val(:∂y))
        ∂p∂z = get_basis_function(ap.bf,Δx,Val(:∂z))
        w, ∂w∂x, ∂w∂y, ∂w∂z = get_kernel_function(ap.kf,Δx,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
        𝝭[i] = p₀ᵀ𝗠⁻¹*p*w
        ∂𝝭∂x[i] = p₀ᵀ∂𝗠⁻¹∂x*p*w + p₀ᵀ𝗠⁻¹*∂p∂x*w + p₀ᵀ𝗠⁻¹*p*∂w∂x
        ∂𝝭∂y[i] = p₀ᵀ∂𝗠⁻¹∂y*p*w + p₀ᵀ𝗠⁻¹*∂p∂y*w + p₀ᵀ𝗠⁻¹*p*∂w∂y
        ∂𝝭∂z[i] = p₀ᵀ∂𝗠⁻¹∂z*p*w + p₀ᵀ𝗠⁻¹*∂p∂z*w + p₀ᵀ𝗠⁻¹*p*∂w∂z
    end
    return 𝝭, ∂𝝭∂x, ∂𝝭∂y, ∂𝝭∂z
end

function cal_moment_matrix!(ap::ReproducingKernel,x::AbstractVector,::Val{:∂1})
    n = get_number_of_basis_function(ap)
    𝗠 = get_moment_matrix(ap,:∂1)
    fill!(𝗠,0.)
    for i in 1:get_number_of_indices(ap)
        xᵢ = get_local_node(ap,i)
        Δx = x - xᵢ
        p = get_basis_function(ap.bf,Δx,Val(:∂1))
        w = get_kernel_function(ap.kf,Δx,Val(:∂1))
        for I in 1:n
            for J in I:n
                𝗠[I,J] += w*p[I]*p[J]
            end
        end
    end
    cholesky!(𝗠)
    U⁻¹ = inverse!(𝗠)
    𝗠⁻¹ = UUᵀ!(U⁻¹)
    return 𝗠⁻¹
end

function cal_moment_matrix!(ap::ReproducingKernel,x::AbstractVector,::Val{:∂1},::Val{:∂x})
    n = get_number_of_basis_function(ap)
    # 𝗠, ∂𝗠∂x = get_moment_matrix(ap,:∂1,:∂x)
    𝗠 = get_moment_matrix(ap,:∂1)
    ∂𝗠∂x = get_moment_matrix(ap,:∂x)
    fill!(𝗠,0.)
    fill!(∂𝗠∂x,0.)
    for i in 1:get_number_of_indices(ap)
        xᵢ = get_local_node(ap,i)
        Δx = x - xᵢ
        # p, ∂p∂x = get_basis_function(ap,Δx,Val(:∂1),Val(:∂x))
        p = get_basis_function(ap.bf,Δx,Val(:∂1))
        ∂p∂x = get_basis_function(ap.bf,Δx,Val(:∂x))
        w, ∂w∂x = get_kernel_function(ap.kf,Δx,Val(:∂1),Val(:∂x))
        for I in 1:n
            for J in I:n
                𝗠[I,J] += w*p[I]*p[J]
                ∂𝗠∂x[I,J] += ∂w∂x*p[I]*p[J] + w*∂p∂x[I]*p[J] + w*p[I]*∂p∂x[J]
            end
        end
    end
    cholesky!(𝗠)
    U⁻¹ = inverse!(𝗠)
    ∂𝗠⁻¹∂x = - UUᵀAUUᵀ!(∂𝗠∂x,U⁻¹)
    𝗠⁻¹ = UUᵀ!(U⁻¹)
    return 𝗠⁻¹, ∂𝗠⁻¹∂x
end

function cal_moment_matrix!(ap::ReproducingKernel,x::AbstractVector,::Val{:∂1},::Val{:∂x},::Val{:∂y},::Val{:∂z})
    n = get_number_of_basis_function(ap)
    # 𝗠, ∂𝗠∂x, ∂𝗠∂y, ∂𝗠∂z = get_moment_matrix(ap,:∂1,:∂x,:∂y,:∂z)
    𝗠 = get_moment_matrix(ap,:∂1)
    ∂𝗠∂x = get_moment_matrix(ap,:∂x)
    ∂𝗠∂y = get_moment_matrix(ap,:∂y)
    ∂𝗠∂z = get_moment_matrix(ap,:∂z)
    fill!(𝗠,0.)
    fill!(∂𝗠∂x,0.)
    fill!(∂𝗠∂y,0.)
    fill!(∂𝗠∂z,0.)
    for i in 1:get_number_of_indices(ap)
        xᵢ = get_local_node(ap,i)
        Δx = x - xᵢ
        # p, ∂p∂x, ∂p∂y, ∂p∂z = get_basis_function(ap,Δx,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
        p = get_basis_function(ap.bf,Δx,Val(:∂1))
        ∂p∂x = get_basis_function(ap.bf,Δx,Val(:∂x))
        ∂p∂y = get_basis_function(ap.bf,Δx,Val(:∂y))
        ∂p∂z = get_basis_function(ap.bf,Δx,Val(:∂z))
        w, ∂w∂x, ∂w∂y, ∂w∂z = get_kernel_function(ap.kf,Δx,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
        for I in 1:n
            for J in I:n
                𝗠[I,J] += w*p[I]*p[J]
                ∂𝗠∂x[I,J] += ∂w∂x*p[I]*p[J] + w*∂p∂x[I]*p[J] + w*p[I]*∂p∂x[J]
                ∂𝗠∂y[I,J] += ∂w∂y*p[I]*p[J] + w*∂p∂y[I]*p[J] + w*p[I]*∂p∂y[J]
                ∂𝗠∂z[I,J] += ∂w∂z*p[I]*p[J] + w*∂p∂z[I]*p[J] + w*p[I]*∂p∂z[J]
            end
        end
    end
    cholesky!(𝗠)
    U⁻¹ = inverse!(𝗠)
    ∂𝗠⁻¹∂x = - UUᵀAUUᵀ!(∂𝗠∂x,U⁻¹)
    ∂𝗠⁻¹∂y = - UUᵀAUUᵀ!(∂𝗠∂y,U⁻¹)
    ∂𝗠⁻¹∂z = - UUᵀAUUᵀ!(∂𝗠∂z,U⁻¹)
    𝗠⁻¹ = UUᵀ!(U⁻¹)
    return 𝗠⁻¹, ∂𝗠⁻¹∂x, ∂𝗠⁻¹∂y, ∂𝗠⁻¹∂z
end

## general functions
# @inline get_basis_function(ap::ReproducingKernel,x::AbstractVector,g::Val) = get_basis_function(ap.bf,x,g)
# @inline get_basis_function(ap::ReproducingKernel,x::AbstractVector,gs::Val...) = (get_basis_function(ap.bf,x,g) for g in gs)
# @inline get_kernel_function(ap::ReproducingKernel,x::AbstractVector,g::Val) = get_kernel_function(ap.kf,x,g)
# @inline get_kernel_function(ap::ReproducingKernel,x::AbstractVector,gs::Val...) = get_kernel_function(ap.kf,x,gs...)
# @inline get_kernel_function(ap::ReproducingKernel,x::AbstractVector,::Val{:∂1}) = get_kernel_function(ap.kf,x,Val(:∂1))
# @inline get_kernel_function(ap::ReproducingKernel,x::AbstractVector,::Val{:∂1},::Val{:∂x}) = get_kernel_function(ap.kf,x,Val(:∂1),Val(:∂x))
# @inline get_kernel_function(ap::ReproducingKernel,x::AbstractVector,::Val{:∂1},::Val{:∂x},::Val{:∂y}) = get_kernel_function(ap.kf,x,Val(:∂1),Val(:∂x),Val(:∂y))
# @inline get_kernel_function(ap::ReproducingKernel,x::AbstractVector,::Val{:∂1},::Val{:∂x},::Val{:∂y},::Val{:∂z}) = get_kernel_function(ap.kf,x,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
# @inline get_kernel(s::Val,r::Float64,gs::Val...) = (get_kernel(s,r,g) for g in gs)
@inline get_moment_matrix(ap::ReproducingKernel,g::Symbol) = ap.bf.𝗠[g]
# @inline get_moment_matrix(ap::ReproducingKernel,gs::Symbol...) = (ap.bf.𝗠[g] for g in gs)
@inline get_shape_function(ap::ReproducingKernel,g::Symbol) = ap.kf.𝝭[g]
# @inline get_shape_function(ap::ReproducingKernel,gs::Symbol...) = (ap.kf.𝝭[g] for g in gs)
@inline get_number_of_basis_function(ap::ReproducingKernel) = ap.bf.𝗠[:∂1].n
@inline get_number_of_shape_functions(ap::ReproducingKernel) = length(ap.kf.𝝭[:∂1])

## spatial partition
struct RegularGrid<:SpatialPartition
    xmin::Vector{Float64}
    dx::Vector{Float64}
    nx::Vector{Int}
    cells::Vector{Set{Int}}
end

# constructions of RegularGrid
function RegularGrid(x::Vector{Node};n::Int=1,γ::Int=1)
    n *= γ
    nₚ  = length(x)
    xmin, xmax = extrema(x[i].x for i in 1:nₚ)
    ymin, ymax = extrema(x[i].y for i in 1:nₚ)
    zmin, zmax = extrema(x[i].z for i in 1:nₚ)
    dx = xmax - xmin
    dy = ymax - ymin
    dz = zmax - zmin
    nd = 0
    pd = 1
    dx > eps() ? (nd += 1;pd *= dx) : dx = 1e-14
    dy > eps() ? (nd += 1;pd *= dy) : dy = 1e-14
    dz > eps() ? (nd += 1;pd *= dz) : dz = 1e-14
    para = (γ*nₚ/pd)^(1/nd)
    nx = ceil(Int, dx * para)
    ny = ceil(Int, dy * para)
    nz = ceil(Int, dz * para)

    cells = Vector{Set{Int}}(undef,nx*ny*nz)
    for i in 1:nx*ny*nz
        cells[i] = Set{Int}()
    end
    for i in 1:nₚ
        ix = floor(Int, (x[i].x - xmin)/dx * nx)
        iy = floor(Int, (x[i].y - ymin)/dy * ny)
        iz = floor(Int, (x[i].z - zmin)/dz * nz)

        ix > nx-1 ? ix = nx-1 : nothing
        iy > ny-1 ? iy = ny-1 : nothing
        iz > nz-1 ? iz = nz-1 : nothing
        for ii in -n:n
            for jj in -n:n
                for kk in -n:n
                    iix = ix + ii
                    iiy = iy + jj
                    iiz = iz + kk

                    iix < 0 ? iix = 0 : nothing
                    iiy < 0 ? iiy = 0 : nothing
                    iiz < 0 ? iiz = 0 : nothing
                    iix > nx-1 ? iix = nx-1 : nothing
                    iiy > ny-1 ? iiy = ny-1 : nothing
                    iiz > nz-1 ? iiz = nz-1 : nothing

                    push!(cells[nx*ny*iiz + nx*iiy + iix + 1], i)
                end
            end
        end
    end
    return RegularGrid([xmin,ymin,zmin],[dx,dy,dz],Int[nx,ny,nz],cells)
end

# actions of RegularGrid
function (rg::RegularGrid)(x::Node)
    ix = floor(Int, (x.x - rg.xmin[1])/rg.dx[1] * rg.nx[1])
    iy = floor(Int, (x.y - rg.xmin[2])/rg.dx[2] * rg.nx[2])
    iz = floor(Int, (x.z - rg.xmin[3])/rg.dx[3] * rg.nx[3])

    ix > rg.nx[1]-1 ? ix = rg.nx[1]-1 : nothing
    iy > rg.nx[2]-1 ? iy = rg.nx[2]-1 : nothing
    iz > rg.nx[3]-1 ? iz = rg.nx[3]-1 : nothing
    return rg.cells[rg.nx[1]*rg.nx[2]*iz + rg.nx[1]*iy + ix + 1]
end

function (rg::RegularGrid)(xs::Node...)
    indices = Set{Int}()
    for x in xs
        union!(indices,rg(x))
    end
    return indices
end
(rg::RegularGrid)(xs::Vector{Node}) = rg(xs...)
