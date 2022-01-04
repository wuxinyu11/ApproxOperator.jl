
## ParametricNode
function get_integration_points(qt::Symbol,::Val{:MFPoint})
    qs = QuadratureRule[qt]
    𝓖 = [MFPoint(q...) for q in qs]
end
# --------------- MFPoint ----------------
struct MFPoint <: ParametricNode
    coordinates::ParametricCoordinates
    w::Float64
    𝝭::Vector{Float64}
    ∂𝝭∂x::Union{Vector{Float64},Nothing}
    ∂𝝭∂y::Union{Vector{Float64},Nothing}
    ∂𝝭∂z::Union{Vector{Float64},Nothing}
    ∂²𝝭∂x²::Union{Vector{Float64},Nothing}
    ∂²𝝭∂x∂y::Union{Vector{Float64},Nothing}
    ∂²𝝭∂y²::Union{Vector{Float64},Nothing}
end

function MFPoint(ξ::T,n::Int,g::Val) where T<:ParametricNode
    𝝭 = zeros(n)
    ∂𝝭∂x = nothing
    ∂𝝭∂y = nothing
    ∂𝝭∂z = nothing
    ∂²𝝭∂x² = nothing
    ∂²𝝭∂x∂y = nothing
    ∂²𝝭∂y² = nothing
    return MFPoint(ξ.coordinates,ξ.w,𝝭,∂𝝭∂x,∂𝝭∂y,∂𝝭∂z,∂²𝝭∂x²,∂²𝝭∂x∂y,∂²𝝭∂y²)
end
function MFPoint(ξ::T,n::Int,gs::Val...) where T<:ParametricNode
    𝝭 = zeros(n)
    ∂𝝭∂x = nothing
    ∂𝝭∂y = nothing
    ∂𝝭∂z = nothing
    ∂²𝝭∂x² = nothing
    ∂²𝝭∂x∂y = nothing
    ∂²𝝭∂y² = nothing
    for g in gs
        if isa(g,Val{:∂x})
            ∂𝝭∂x = zeros(n)
        elseif isa(g,Val{:∂y})
            ∂𝝭∂y = zeros(n)
        elseif isa(g,Val{:∂z})
            ∂𝝭∂z = zeros(n)
        elseif isa(g,Val{:∂x²})
            ∂²𝝭∂x² = zeros(n)
        elseif isa(g,Val{:∂x∂y})
            ∂²𝝭∂x∂y = zeros(n)
        elseif isa(g,Val{:∂y²})
            ∂²𝝭∂y² = zeros(n)
        end
    end

    return MFPoint(ξ.coordinates,ξ.w,𝝭,∂𝝭∂x,∂𝝭∂y,∂𝝭∂z,∂²𝝭∂x²,∂²𝝭∂x∂y,∂²𝝭∂y²)
end

# action
get_shape_functions(::Approximator,ξ::MFPoint,::Val{:∂1}) = ξ.𝝭
get_shape_functions(::Approximator,ξ::MFPoint,::Val{:∂x}) = ξ.∂𝝭∂x
get_shape_functions(::Approximator,ξ::MFPoint,::Val{:∂y}) = ξ.∂𝝭∂y
get_shape_functions(::Approximator,ξ::MFPoint,::Val{:∂z}) = ξ.∂𝝭∂z
get_shape_functions(::Approximator,ξ::MFPoint,::Val{:∂x²}) = ξ.∂²𝝭∂x²
get_shape_functions(::Approximator,ξ::MFPoint,::Val{:∂y²}) = ξ.∂²𝝭∂y²

# --------------- RKPoint ----------------
struct RKPoint <: ParametricNode
    coordinates::ParametricCoordinates
    wⁱ::Float64
    wᵇ::Float64
    𝝭::Vector{Float64}
end
