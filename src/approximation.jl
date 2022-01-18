
## AbstractPoi
@inline getx(ap::AbstractPoi,::Node) = (ap.𝓒[1].x,ap.𝓒[1].y,ap.𝓒[1].z)
@inline getw(ap::AbstractPoi,::Node) = 1.0
# -------------- Poi1 --------------
struct Poi1<:AbstractPoi
    𝓒::Vector{Node}
    𝓖::Vector{Node}
end
get𝝭(::Poi1,::Node) = 1.0

## AbstractSeg
@inline getx(ap::A,ξ::Node) where A<:AbstractSeg = getx(ap,ξ.ξ)
@inline function getx(ap::A,ξ::Float64) where A<:AbstractSeg
    x₁ = ap.𝓒[1].x
    y₁ = ap.𝓒[1].y
    z₁ = ap.𝓒[1].z
    x₂ = ap.𝓒[2].x
    y₂ = ap.𝓒[2].y
    z₂ = ap.𝓒[2].z
    N₁ = 0.5*(1-ξ)
    N₂ = 0.5*(1+ξ)
    return (x₁*N₁+x₂*N₂,y₁*N₁+y₂*N₂,z₁*N₁+z₂*N₂)
end

@inline getw(ap::A,ξ::Node) where A<:AbstractSeg = 0.5*ap.L*ξ.w

# ---------------- Seg2 -------------------
struct Seg2<:AbstractSeg
    𝓒::Vector{Node}
    𝓖::Vector{Node}
    L::Float64
end

# constructions of Seg2
function Seg2(𝓒::Vector{Node},𝓖::Vector{Node})
    x₁ = 𝓒[1].x
    y₁ = 𝓒[1].y
    x₂ = 𝓒[2].x
    y₂ = 𝓒[2].y
    L = ((x₂-x₁)^2+(y₂-y₁)^2)^0.5
    return Seg2(𝓒,𝓖,L)
end

# actions for Seg2
@inline get𝝭(ap::Seg2,ξ::Node) = get𝝭(ap,ξ.ξ)
@inline get𝝭(ap::Seg2,ξ::Float64) = (0.5*(1-ξ),0.5*(1+ξ))
@inline get∂𝝭∂x(ap::Seg2,::Node) = (-1.0/ap.L,1.0/ap.L)
@inline get∂𝝭∂x(ap::Seg2,::Float64) = (-1.0/ap.L,1.0/ap.L)
@inline get∂𝝭∂y(ap::Seg2,::Node) = (0.0,0.0)
@inline get∂𝝭∂z(ap::Seg2,::Node) = (0.0,0.0)
@inline get∇𝝭(ap::Seg2,ξ::Node) = (get𝝭(ap,ξ),get∂𝝭∂x(ap,ξ),(0.0,0.0),(0.0,0.0))

##
struct Tri3
    fields
end


##
struct Quad
    fields
end

## PoiN
struct PoiN{T,𝒑,𝑠,𝜙}<:ReproducingKernel{T,𝒑,𝑠,𝜙}
    𝓒::Vector{Node}
    𝓖::Vector{T}
    𝗠::Dict{Symbol,SymMat}
    𝝭::Dict{Symbol,Vector{Float64}}
    type::Tuple{Val{𝒑},Val{𝑠},Val{𝜙}}
end

PoiN(𝓒::Vector{Node},𝓖::Vector{T},𝗠::Dict{Symbol,SymMat},𝝭::Dict{Symbol,Vector{Float64}},𝒑::Symbol,𝑠::Symbol,𝜙::Symbol) where T<:AbstractNode = PoiN(𝓒,𝓖,𝗠,𝝭,(Val(𝒑),Val(𝑠),Val(𝜙)))

@inline getx(ap::PoiN,::AbstractNode) = (ap.𝓒[1].x,ap.𝓒[1].y,ap.𝓒[1].z)
@inline getw(ap::PoiN,::Node) = 1.0

## SegN
struct SegN{T,𝒑,𝑠,𝜙}<:ReproducingKernel{T,𝒑,𝑠,𝜙}
    𝓒::Vector{Node}
    𝓖::Vector{T}
    𝗠::Dict{Symbol,SymMat}
    𝝭::Dict{Symbol,Vector{Float64}}
    type::Tuple{Val{𝒑},Val{𝑠},Val{𝜙}}
    L::Float64
end

function SegN(𝓒::Vector{Node},𝓖::Vector{T},𝗠::Dict{Symbol,SymMat},𝝭::Dict{Symbol,Vector{Float64}},𝒑::Symbol,𝑠::Symbol,𝜙::Symbol) where T<:AbstractNode
    x₁ = 𝓒[1].x
    y₁ = 𝓒[1].y
    x₂ = 𝓒[2].x
    y₂ = 𝓒[2].y
    L = ((x₂-x₁)^2+(y₂-y₁)^2)^0.5

    return SegN(𝓒,𝓖,𝗠,𝝭,(Val(𝒑),Val(𝑠),Val(𝜙)),L)
end

@inline getx(ap::SegN,ξ::T) where T<:AbstractNode = getx(ap,ξ.ξ)
@inline function getx(ap::SegN,ξ::Float64)
    x₁ = ap.𝓒[1].x
    y₁ = ap.𝓒[1].y
    z₁ = ap.𝓒[1].z
    x₂ = ap.𝓒[2].x
    y₂ = ap.𝓒[2].y
    z₂ = ap.𝓒[2].z
    N₁ = 0.5*(1-ξ)
    N₂ = 0.5*(1+ξ)
    return (x₁*N₁+x₂*N₂,y₁*N₁+y₂*N₂,z₁*N₁+z₂*N₂)
end

@inline getw(ap::SegN,ξ::T) where T<:AbstractNode = 0.5*ap.L*ξ.w

##
struct TriN
    fields
end

struct QuadN
    fields
end

function get𝝭(ap::T,ξ::Node) where T<:ReproducingKernel
    𝓒 = ap.𝓒
    𝝭 = ap.𝝭[:∂1]
    x = getx(ap,ξ)
    𝒑₀ᵀ𝗠⁻¹ = cal𝗠!(ap,x)
    for i in 1:length(𝓒)
        xᵢ = 𝓒[i]
        Δx = x - xᵢ
        𝒑 = get𝒑(ap,Δx)
        𝜙 = get𝜙(ap,xᵢ,Δx)
        𝝭[i] = 𝒑₀ᵀ𝗠⁻¹*𝒑*𝜙
    end
    return 𝝭
end

function get∂𝝭∂x(ap::T,ξ::Node) where T<:ReproducingKernel
    𝓒 = ap.𝓒
    𝝭 = ap.𝝭[:∂1]
    ∂𝝭∂x = ap.𝝭[:∂x]
    x = getx(ap,ξ)
    p₀ᵀ𝗠⁻¹, p₀ᵀ∂𝗠⁻¹∂x = cal∂𝗠∂x!(ap,x)
    for i in 1:length(𝓒)
        xᵢ = 𝓒[i]
        Δx = x - xᵢ
        𝒑, ∂𝒑∂x = get∂𝒑∂x(ap,Δx)
        𝜙, ∂𝜙∂x = get∂𝜙∂x(ap,xᵢ,Δx)
        𝝭[i] = 𝒑₀ᵀ𝗠⁻¹*𝒑*𝜙
        ∂𝝭∂x[i] = 𝒑₀ᵀ∂𝗠⁻¹∂x*𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂x*𝜙 + 𝒑₀ᵀ𝗠⁻¹*𝒑*∂𝜙∂x
    end
    return 𝝭, ∂𝝭∂x
end

function get∇𝝭(ap::T,ξ::Node) where T<:ReproducingKernel
    𝓒 = ap.𝓒
    𝝭 = ap.𝝭[:∂1]
    ∂𝝭∂x = ap.𝝭[:∂x]
    ∂𝝭∂y = ap.𝝭[:∂y]
    ∂𝝭∂z = ap.𝝭[:∂z]
    x = getx(ap,ξ)
    𝒑₀ᵀ𝗠⁻¹, 𝒑₀ᵀ∂𝗠⁻¹∂x, 𝒑₀ᵀ∂𝗠⁻¹∂y, 𝒑₀ᵀ∂𝗠⁻¹∂z= cal∇𝗠!(ap,x)
    for i in 1:length(𝓒)
        xᵢ = 𝓒[i]
        Δx = x - xᵢ
        𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂𝒑∂z = get∇𝒑(ap,Δx)
        𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂𝜙∂z = get∇𝜙(ap,xᵢ,Δx)
        𝝭[i] = 𝒑₀ᵀ𝗠⁻¹*𝒑*𝜙
        ∂𝝭∂x[i] = 𝒑₀ᵀ∂𝗠⁻¹∂x*𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂x*𝜙 + 𝒑₀ᵀ𝗠⁻¹*𝒑*∂𝜙∂x
        ∂𝝭∂y[i] = 𝒑₀ᵀ∂𝗠⁻¹∂y*𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂y*𝜙 + 𝒑₀ᵀ𝗠⁻¹*𝒑*∂𝜙∂y
        ∂𝝭∂z[i] = 𝒑₀ᵀ∂𝗠⁻¹∂z*𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂z*𝜙 + 𝒑₀ᵀ𝗠⁻¹*𝒑*∂𝜙∂z
    end
    return 𝝭, ∂𝝭∂x, ∂𝝭∂y, ∂𝝭∂z
end
