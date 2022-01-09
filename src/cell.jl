
##
@inline getproperty(ap::T,f::Symbol) where T<:Approximator = getdata(ap,Val(f))
@inline getdata(ap::T,::Val{:𝓒}) where T<:Approximator = getfield(ap,:𝓒)
@inline getdata(ap::T,::Val{:𝓖}) where T<:Approximator = getfield(ap,:𝓖)
@inline getdata(ap::T,::Val{:𝝭}) where T<:Approximator = (ξ::Point)->ap.𝝭(ξ.coordinates)
@inline function getdata(ap::T,::Val{:coordinates}) where T<:Approximator
    return (ξ)->(sum(ap.𝝭(ξ)[i]*ap.𝓒[i].x for i in 1:length(ap.𝓒)),sum(ap.𝝭(ξ)[j]*ap.𝓒[j].y for j in 1:length(ap.𝓒)),sum(ap.𝝭(ξ)[k]*ap.𝓒[k].z for k in 1:length(ap.𝓒)))
end

## AbstractPoi
@inline getdata(ap::T,::Val{:J}) where T<:AbstractPoi = (::Any)->1.0

# ---------------- Poi1 ----------------
struct Poi1{T}<:AbstractPoi where T<:ParametricNode
    𝓒::Node
    𝓖::Vector{T}
end
@inline getdata(::Poi1,::Val{:𝝭}) = (::Any)->1.0

## AbstractSeg
@inline getdata(ap::T,::Val{:L}) where T<:AbstractSeg = getfield(ap,:L)
@inline getdata(ap::T,::Val{:J}) where T<:AbstractSeg = (::Any)->0.5*ap.L

# ---------------- Seg2 -------------------
struct Seg2{T}<:AbstractSeg where T<:ParametricNode
    𝓒::Vector{Node}
    𝓖::Vector{T}
    L::Float64
end

# constructions of Seg2
function Seg2(𝓒::Vector{Node},𝓖::Vector{T}) where T<:ParametricNode
    x₁ = 𝓒[1].x
    y₁ = 𝓒[1].y
    x₂ = 𝓒[2].x
    y₂ = 𝓒[2].y
    L = ((x₂-x₁)^2+(y₂-y₁)^2)^0.5
    return Seg2(𝓒,𝓖,L)
end

# actions for Seg2
@inline function getdata(ap::Seg2,::Val{:𝝭})
    @inline get𝝭(ξ::Float64) = ((1.0-ξ)*0.5,(1.0+ξ)*0.5)
    @inline get𝝭(ξ::Point) = ((1.0-ξ.ξ)*0.5,(1.0+ξ.ξ)*0.5)
    return get𝝭
end
@inline getdata(ap::Seg2,::Val{:∂𝝭∂x}) = (::Any)->(-1.0/ap.L,1.0/ap.L)
@inline getdata(  ::Seg2,::Val{:∂𝝭∂y}) = (::Any)->(0.,0.)
@inline getdata(  ::Seg2,::Val{:∂𝝭∂z}) = (::Any)->(0.,0.)
@inline getdata(ap::Seg2,::Val{:∇𝝭}) = (ξ)->(ap.∂𝝭∂x(ξ),ap.∂𝝭∂y(ξ),ap.∂𝝭∂z(ξ))
