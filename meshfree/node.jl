
## ParametricNode
# --------------- RKPoint ----------------
struct RKPoint <: ParametricNode
    coordinates::ParametricCoordinates
    wⁱ::Float64
    wᵇ::Float64
    𝝭::Vector{Float64}
end

# --------------- MFPoint ---------------
struct MFPoint <: ParametricNode
    coordinates::ParametricCoordinates
    w::Float64
    bf::BasisFunction
    kf::KernelFunction
end
