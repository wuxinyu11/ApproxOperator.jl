
## ParametricNode
struct MFPoint <: ParametricNode
    coordinates::ParametricCoordinates
    𝝭::Vector{Float64}
    ∂𝝭∂x::Union{Vector{Float64},Nothing}
    ∂𝝭∂y::Union{Vector{Float64},Nothing}
    ∂𝝭∂z::Union{Vector{Float64},Nothing}
    ∂²𝝭∂x²::Union{Vector{Float64},Nothing}
    ∂²𝝭∂xy::Union{Vector{Float64},Nothing}
    ∂²𝝭∂y²::Union{Vector{Float64},Nothing}
end

# action
get_shape_functions(::Approximator,ξ::MFPoint,::Val{:∂1}) = ξ.𝝭
get_shape_functions(::Approximator,ξ::MFPoint,::Val{:∂x}) = ξ.∂𝝭∂x
get_shape_functions(::Approximator,ξ::MFPoint,::Val{:∂y}) = ξ.∂𝝭∂y
get_shape_functions(::Approximator,ξ::MFPoint,::Val{:∂z}) = ξ.∂𝝭∂z
get_shape_functions(::Approximator,ξ::MFPoint,::Val{:∂x²}) = ξ.∂²𝝭∂x²
get_shape_functions(::Approximator,ξ::MFPoint,::Val{:∂y²}) = ξ.∂²𝝭∂y²
