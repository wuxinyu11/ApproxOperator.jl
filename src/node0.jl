ParametricCoordinates = Union{Float64,NTuple{2,Float64},NTuple{3,Float64}}
## PhysicalNode
@inline +(n::NTuple{3,Float64},m::NTuple{3,Float64}) = (n[1]+m[1], n[2]+m[2], n[3]+m[3])
@inline -(n::NTuple{3,Float64},m::NTuple{3,Float64}) = (n[1]-m[1], n[2]-m[2], n[3]-m[3])
@inline *(c::Float64,n::NTuple{3,Float64}) = (c*n[1], c*n[2], c*n[3])
# ------------- Node -------------
struct Node <: PhysicalNode
    coordinates::NTuple{3,Float64}
end
Node(x::Float64,y::Float64,z::Float64) = Node((x,y,z))

## ParametricNode
# --------------- Gauss integration point ----------------
struct GaussPoint <: ParametricNode
    coordinates::ParametricCoordinates
    w::Float64
end
GaussPoint(ξ₁::Float64,ξ₂::Float64,w::Float64) = GaussPoint((ξ₁,ξ₂),w)
GaussPoint(ξ₁::Float64,ξ₂::Float64,ξ₃::Float64,w::Float64) = GaussPoint((ξ₁,ξ₂,ξ₃),w)

# --------------- SSPoint ----------------
struct SSPoint <: ParametricNode
    coordinates::ParametricCoordinates
    w::Float64
    s::Int
    𝝭::SparseShapePool
end

function SSPoint(ξ::T,s::Int,𝝭::SparseShapePool) where T<:ParametricNode
    return SSPoint(ξ.coordinates,ξ.w,s,𝝭)
end

# action
get_shape_functions(::T,ξ::SSPoint,::Val{:∂1}) where T<:Approximator = SparseShape(ξ.𝝭,ξ.s,Val(:∂1))
get_shape_functions(::T,ξ::SSPoint,::Val{:∂x}) where T<:Approximator = SparseShape(ξ.𝝭,ξ.s,Val(:∂x))
get_shape_functions(::T,ξ::SSPoint,::Val{:∂y}) where T<:Approximator = SparseShape(ξ.𝝭,ξ.s,Val(:∂y))
get_shape_functions(::T,ξ::SSPoint,::Val{:∂z}) where T<:Approximator = SparseShape(ξ.𝝭,ξ.s,Val(:∂z))
get_shape_functions(::T,ξ::SSPoint,::Val{:∂x²}) where T<:Approximator = SparseShape(ξ.𝝭,ξ.s,Val(:∂x²))
get_shape_functions(::T,ξ::SSPoint,::Val{:∂x∂y}) where T<:Approximator = SparseShape(ξ.𝝭,ξ.s,Val(:∂x∂y))
get_shape_functions(::T,ξ::SSPoint,::Val{:∂y²}) where T<:Approximator = SparseShape(ξ.𝝭,ξ.s,Val(:∂y²))

## Sparse shape function storge
struct SparseShapePool
    index::Vector{Int}
    𝝭::Vector{Float64}
    ∂𝝭∂x::Union{Vector{Float64},Nothing}
    ∂𝝭∂y::Union{Vector{Float64},Nothing}
    ∂𝝭∂z::Union{Vector{Float64},Nothing}
    ∂²𝝭∂x²::Union{Vector{Float64},Nothing}
    ∂²𝝭∂x∂y::Union{Vector{Float64},Nothing}
    ∂²𝝭∂y²::Union{Vector{Float64},Nothing}
end

function SparseShapePool(n::Int,gs::Val...)
    𝝭 = nothing
    ∂𝝭∂x = nothing
    ∂𝝭∂y = nothing
    ∂𝝭∂z = nothing
    ∂²𝝭∂x² = nothing
    ∂²𝝭∂x∂y = nothing
    ∂²𝝭∂y² = nothing
    for g in gs
        if isa(g,Val{:∂1})
            𝝭 = Vector{Float64}[]
        elseif isa(g,Val{:∂x})
            ∂𝝭∂x = Vector{Float64}[]
        elseif isa(g,Val{:∂y})
            ∂𝝭∂y = Vector{Float64}[]
        elseif isa(g,Val{:∂z})
            ∂𝝭∂z = Vector{Float64}[]
        elseif isa(g,Val{:∂x²})
            ∂²𝝭∂x² = Vector{Float64}[]
        elseif isa(g,Val{:∂x∂y})
            ∂²𝝭∂x∂y = Vector{Float64}[]
        elseif isa(g,Val{:∂y²})
            ∂²𝝭∂y² = Vector{Float64}[]
        end
    end
    return SparseShapePool(zeros(n+1),𝝭,∂𝝭∂x,∂𝝭∂y,∂𝝭∂z,∂²𝝭∂x²,∂²𝝭∂x∂y,∂²𝝭∂y²)
end

(sp::SparseShapePool)(n::Int,::Val{:∂1}) = sp.𝝭[n]
(sp::SparseShapePool)(n::Int,::Val{:∂x}) = sp.∂𝝭∂x[n]
(sp::SparseShapePool)(n::Int,::Val{:∂y}) = sp.∂𝝭∂y[n]
(sp::SparseShapePool)(n::Int,::Val{:∂z}) = sp.∂𝝭∂z[n]
(sp::SparseShapePool)(n::Int,::Val{:∂x²}) = sp.∂²𝝭∂x²[n]
(sp::SparseShapePool)(n::Int,::Val{:∂x∂y}) = sp.∂²𝝭∂x∂y[n]
(sp::SparseShapePool)(n::Int,::Val{:∂y²}) = sp.∂²𝝭∂y²[n]

struct SparseShape
    p::SparseShapePool
    n::Int
    t::Val
end

function getindex(ss::SparseShape,i::Int)
    n = ss.p.index[ss.n]+i
    return ss.p(n,ss.t)
end
## Actions
getindex(x::T,i::Int) where T<:AbstractNode = x.coordinates[i]
