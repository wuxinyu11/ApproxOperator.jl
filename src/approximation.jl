function set𝒏(aps::Vector{A}) where A<:Approximator
    for ap in aps
        set𝒏(ap)
    end
end
function set𝒏(ap::A) where A<:AbstractSeg
    𝓖 = ap.𝓖
    for ξ in 𝓖
        ξ.n₁ = get𝒏(ap,ξ)
    end
end
function set𝒏(ap::A) where A<:AbstractTri
    𝓖 = ap.𝓖
    for ξ in 𝓖
        ξ.n₁, ξ.n₂ = get𝒏(ap,ξ)
    end
end
## AbstractPoi
@inline getx(ap::AbstractPoi,::AbstractNode) = (ap.𝓒[1].x,ap.𝓒[1].y,ap.𝓒[1].z)
@inline getw(ap::AbstractPoi,::AbstractNode) = 1.0
# -------------- Poi1 --------------
struct Poi1<:AbstractPoi
    𝓒::Vector{Node}
    𝓖::Vector{Node}
end
function Poi1(i::Int,data::Dict{Symbol,Vector{Float64}}) where T<:AbstractNode
    𝓒 = [Node(i,data)]
    𝓖 = Node[]
    return Poi1(𝓒,𝓖)
end
get𝝭(::Poi1,::Node) = 1.0

## AbstractSeg
@inline getx(ap::A,ξ::T) where {A<:AbstractSeg,T<:AbstractNode} = getx(ap,ξ.ξ)
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

@inline getw(ap::A,ξ::T) where {A<:AbstractSeg,T<:AbstractNode} = 0.5*ap.L*ξ.w
@inline get𝒏(ap::A,ξ::T) where {A<:AbstractSeg,T<:AbstractNode} = get𝒏(ap,ξ.ξ)
function get𝒏(ap::A,ξ::Float64) where A<:AbstractSeg
    n₁ = 0.0
    n₁ += ξ == -1.0 ?  1.0 : 0.0
    n₁ += ξ ==  1.0 ? -1.0 : 0.0
end

# ---------------- Seg2 -------------------
struct Seg2<:AbstractSeg
    𝓒::Vector{Node}
    𝓖::Vector{Node}
    L::Float64
end
function Seg2(i::Int,j::Int,data::Dict{Symbol,Vector{Float64}})
    𝓒 = [Node(i,data),Node(j,data)]
    𝓖 = Node[]
    return Seg2(𝓒,𝓖)
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

## AbstractTri
@inline getx(ap::A,ξ::T) where {A<:AbstractTri,T<:AbstractNode} = getx(ap,ξ.ξ,ξ.η)
@inline function getx(ap::A,ξ::Float64,η::Float64) where A<:AbstractTri
    x₁ = ap.𝓒[1].x
    y₁ = ap.𝓒[1].y
    z₁ = ap.𝓒[1].z
    x₂ = ap.𝓒[2].x
    y₂ = ap.𝓒[2].y
    z₂ = ap.𝓒[2].z
    x₃ = ap.𝓒[3].x
    y₃ = ap.𝓒[3].y
    z₃ = ap.𝓒[3].z
    N₁ = ξ
    N₂ = η
    N₃ = 1.0-ξ-η
    return (x₁*N₁+x₂*N₂+x₃*N₃,y₁*N₁+y₂*N₂+y₃*N₃,z₁*N₁+z₂*N₂+z₃*N₃)
end
@inline getw(ap::A,ξ::T) where {A<:AbstractTri,T<:AbstractNode} = ap.A*ξ.w
@inline get𝒏(ap::A,ξ::T) where {A<:AbstractTri,T<:AbstractNode} = get𝒏(ap,ξ.ξ,ξ.η)
function get𝒏(ap::A,ξ::Float64,η::Float64) where A<:AbstractTri
    n₁ = 0
    n₂ = 0
    x₁ = ap.𝓒[1].x
    y₁ = ap.𝓒[1].y
    x₂ = ap.𝓒[2].x
    y₂ = ap.𝓒[2].y
    x₃ = ap.𝓒[3].x
    y₃ = ap.𝓒[3].y
    γ = 1.0-ξ-η
    n₁ += ξ == 0.0 ? y₃-y₂ : 0.0
    n₁ += η == 0.0 ? y₁-y₃ : 0.0
    n₁ += γ == 0.0 ? y₂-y₁ : 0.0
    n₂ += ξ == 0.0 ? x₂-x₃ : 0.0
    n₂ += η == 0.0 ? x₃-x₁ : 0.0
    n₂ += γ == 0.0 ? x₁-x₂ : 0.0
    return n₁,n₂
end

struct Tri3<:AbstractTri
    𝓒::Vector{Node}
    𝓖::Vector{Node}
    A::Float64
end

function Tri3(i::Int,j::Int,k::Int,data::Dict{Symbol,Vector{Float64}})
    𝓒 = [Node(i,data),Node(j,data),Node(k,data)]
    𝓖 = Node[]
    return Tri3(𝓒,𝓖)
end

# constructions of Tri3
function Tri3(𝓒::Vector{Node},𝓖::Vector{Node})
    x₁ = 𝓒[1].x
    y₁ = 𝓒[1].y
    z₁ = 𝓒[1].z
    x₂ = 𝓒[2].x
    y₂ = 𝓒[2].y
    z₂ = 𝓒[2].z
    x₃ = 𝓒[3].x
    y₃ = 𝓒[3].y
    z₃ = 𝓒[3].z
    A₁ = 0.5*(y₁*z₂+y₂*z₃+y₃*z₁-y₂*z₁-y₃*z₂-y₁*z₃)
    A₂ = 0.5*(z₁*x₂+z₂*x₃+z₃*x₁-z₂*x₁-z₃*x₂-z₁*x₃)
    A₃ = 0.5*(x₁*y₂+x₂*y₃+x₃*y₁-x₂*y₁-x₃*y₂-x₁*y₃)
    A = (A₁^2 + A₂^2 + A₃^2)^0.5
    return Tri3(𝓒,𝓖,A)
end

## AbstractQuad
function getw(ap::AbstractQuad,ξ::AbstractNode)
    J₁₁,J₂₁,J₁₂,J₂₂ = getJ(ap,ξ)
    return (J₁₁*J₂₂-J₂₁*J₁₂)*ξ.w
end
function getJ(ap::AbstractQuad,ξ::AbstractNode)
    x₁ = ap.𝓒[1].x
    x₂ = ap.𝓒[2].x
    x₃ = ap.𝓒[3].x
    x₄ = ap.𝓒[4].x
    y₁ = ap.𝓒[1].y
    y₂ = ap.𝓒[2].y
    y₃ = ap.𝓒[3].y
    y₄ = ap.𝓒[4].y
    ∂N₁∂ξ,∂N₂∂ξ,∂N₃∂ξ,∂N₄∂ξ = get∂𝝭∂ξ(ap,ξ)
    ∂N₁∂η,∂N₂∂η,∂N₃∂η,∂N₄∂η = get∂𝝭∂η(ap,ξ)
    J₁₁ = ∂N₁∂ξ*x₁ + ∂N₂∂ξ*x₂ + ∂N₃∂ξ*x₃ + ∂N₄∂ξ*x₄
    J₁₂ = ∂N₁∂η*x₁ + ∂N₂∂η*x₂ + ∂N₃∂η*x₃ + ∂N₄∂η*x₄
    J₂₁ = ∂N₁∂ξ*y₁ + ∂N₂∂ξ*y₂ + ∂N₃∂ξ*y₃ + ∂N₄∂ξ*y₄
    J₂₂ = ∂N₁∂η*y₁ + ∂N₂∂η*y₂ + ∂N₃∂η*y₃ + ∂N₄∂η*y₄
    return J₁₁,J₂₁,J₁₂,J₂₂
end
function getx(ap::AbstractQuad,ξ::AbstractNode)
    N₁,N₂,N₃,N₄ = get𝝭(ap,ξ)
    x₁ = ap.𝓒[1].x
    x₂ = ap.𝓒[2].x
    x₃ = ap.𝓒[3].x
    x₄ = ap.𝓒[4].x
    y₁ = ap.𝓒[1].y
    y₂ = ap.𝓒[2].y
    y₃ = ap.𝓒[3].y
    y₄ = ap.𝓒[4].y
    z₁ = ap.𝓒[1].z
    z₂ = ap.𝓒[2].z
    z₃ = ap.𝓒[3].z
    z₄ = ap.𝓒[4].z
    return (N₁*x₁+N₂*x₂+N₃*x₃+N₄*x₄, N₁*y₁+N₂*y₂+N₃*y₃+N₄*y₄, N₁*z₁+N₂*z₂+N₃*z₃+N₄*z₄)
end

struct Quad
    𝓒::Vector{Node}
    𝓖::Vector{Node}
end

function Quad(i::Int,j::Int,k::Int,l::Int,data::Dict{Symbol,Vector{Float64}})
    𝓒 = [Node(i,data),Node(j,data),Node(k,data),Node(l,data)]
    𝓖 = Node[]
    return Quad(𝓒,𝓖)
end

get𝝭(ap::Quad,ξ::AbstractNode) = get𝝭(ap,ξ.ξ,ξ.η)
get∂𝝭∂ξ(ap::Quad,ξ::AbstractNode) = get∂𝝭∂ξ(ap,ξ.η)
get∂𝝭∂η(ap::Quad,ξ::AbstractNode) = get∂𝝭∂η(ap,ξ.ξ)

function get𝝭(ap::Quad,ξ::Float64,η::Float64)
    N₁ = 0.25*(1.0-ξ)*(1.0-η)
    N₂ = 0.25*(1.0+ξ)*(1.0-η)
    N₃ = 0.25*(1.0+ξ)*(1.0+η)
    N₄ = 0.25*(1.0-ξ)*(1.0+η)
    return (N₁,N₂,N₃,N₄)
end
function get∂𝝭∂ξ(ap::Quad,η::Float64)
    ∂N₁∂ξ = - 0.25*(1-η)
    ∂N₂∂ξ =   0.25*(1-η)
    ∂N₃∂ξ =   0.25*(1+η)
    ∂N₄∂ξ = - 0.25*(1+η)
    return (∂N₁∂ξ,∂N₂∂ξ,∂N₃∂ξ,∂N₄∂ξ)
end
function get∂𝝭∂η(ap::Quad,ξ::Float64)
    ∂N₁∂η = - 0.25*(1-ξ)
    ∂N₂∂η = - 0.25*(1+ξ)
    ∂N₃∂η =   0.25*(1+ξ)
    ∂N₄∂η =   0.25*(1-ξ)
    return (∂N₁∂η,∂N₂∂η,∂N₃∂η,∂N₄∂η)
end
function get∂𝝭∂x∂𝝭∂y(ap::Quad,ξ::AbstractNode)
    x₁ = ap.𝓒[1].x
    x₂ = ap.𝓒[2].x
    x₃ = ap.𝓒[3].x
    x₄ = ap.𝓒[4].x
    y₁ = ap.𝓒[1].y
    y₂ = ap.𝓒[2].y
    y₃ = ap.𝓒[3].y
    y₄ = ap.𝓒[4].y
    ∂N₁∂ξ,∂N₂∂ξ,∂N₃∂ξ,∂N₄∂ξ = get∂𝝭∂ξ(ap,ξ)
    ∂N₁∂η,∂N₂∂η,∂N₃∂η,∂N₄∂η = get∂𝝭∂η(ap,ξ)
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
get∇𝝭(ap::Quad,ξ::AbstractNode) = get𝝭(ap,ξ),get∂𝝭∂x∂𝝭∂y(ap,ξ)...,(0.0,0.0,0.0,0.0)

## PoiN
struct PoiN{T,𝒑,𝑠,𝜙}<:ReproducingKernel{T,𝒑,𝑠,𝜙}
    𝓒::Vector{Node}
    𝓖::Vector{T}
    𝗠::Dict{Symbol,SymMat}
    𝝭::Dict{Symbol,Vector{Float64}}
    type::Tuple{Val{𝒑},Val{𝑠},Val{𝜙}}
end

function PoiN{T,𝒑,𝑠,𝜙}(𝓒::Vector{Node},𝗠::Dict{Symbol,SymMat},𝝭::Dict{Symbol,Vector{Float64}}) where {T<:AbstractNode,𝒑,𝑠,𝜙}
    𝓖 = T[]
    return PoiN{T,𝒑,𝑠,𝜙}(𝓒,𝓖,𝗠,𝝭,(Val(𝒑),Val(𝑠),Val(𝜙)))
end
function PoiN{T,𝒑,𝑠,𝜙}(i::Int,data::Dict{Symbol,Vector{Float64}},𝗠::Dict{Symbol,SymMat},𝝭::Dict{Symbol,Vector{Float64}}) where {T<:AbstractNode,𝒑,𝑠,𝜙}
    𝓒 = [Node(i,data)]
    return PoiN{T,𝒑,𝑠,𝜙}(𝓒,𝗠,𝝭)
end

@inline getx(ap::PoiN,::AbstractNode) = (ap.𝓒[1].x,ap.𝓒[1].y,ap.𝓒[1].z)
@inline getw(ap::PoiN,::AbstractNode) = 1.0

## SegN
struct SegN{T,𝒑,𝑠,𝜙}<:ReproducingKernel{T,𝒑,𝑠,𝜙}
    𝓒::Vector{Node}
    𝓖::Vector{T}
    𝗠::Dict{Symbol,SymMat}
    𝝭::Dict{Symbol,Vector{Float64}}
    type::Tuple{Val{𝒑},Val{𝑠},Val{𝜙}}
    L::Float64
end

function SegN{T,𝒑,𝑠,𝜙}(𝓒::Vector{Node},𝗠::Dict{Symbol,SymMat},𝝭::Dict{Symbol,Vector{Float64}}) where {T<:AbstractNode,𝒑,𝑠,𝜙}
    x₁ = 𝓒[1].x
    y₁ = 𝓒[1].y
    x₂ = 𝓒[2].x
    y₂ = 𝓒[2].y
    L = ((x₂-x₁)^2+(y₂-y₁)^2)^0.5
    𝓖 = T[]

    return SegN(𝓒,𝓖,𝗠,𝝭,(Val(𝒑),Val(𝑠),Val(𝜙)),L)
end

function SegN{T,𝒑,𝑠,𝜙}(i::Int,j::Int,data::Dict{Symbol,Vector{Float64}},𝗠::Dict{Symbol,SymMat},𝝭::Dict{Symbol,Vector{Float64}}) where {T<:AbstractNode,𝒑,𝑠,𝜙}
    𝓒 = [Node(i,data),Node(j,data)]
    return SegN{T,𝒑,𝑠,𝜙}(𝓒,𝗠,𝝭)
end

@inline getx(ap::SegN,ξ::AbstractNode) = getx(ap,ξ.ξ)
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
function get𝒏(ap::SegN,ξ::T) where T<:AbstractNode
    n₁ = 0.0
    n₁ += ξ.ξ ==  1.0 ?  1.0 : 0.0
    n₁ += ξ.ξ == -1.0 ? -1.0 : 0.0
end

##
struct TriN{T,𝒑,𝑠,𝜙}<:ReproducingKernel{T,𝒑,𝑠,𝜙}
    𝓒::Vector{Node}
    𝓖::Vector{T}
    𝗠::Dict{Symbol,SymMat}
    𝝭::Dict{Symbol,Vector{Float64}}
    type::Tuple{Val{𝒑},Val{𝑠},Val{𝜙}}
    A::Float64
end

struct TetN{T,𝒑,𝑠,𝜙}<:ReproducingKernel{T,𝒑,𝑠,𝜙}
    𝓒::Vector{Node}
    𝓖::Vector{T}
    𝗠::Dict{Symbol,SymMat}
    𝝭::Dict{Symbol,Vector{Float64}}
    type::Tuple{Val{𝒑},Val{𝑠},Val{𝜙}}
    V::Float64
end

## get shape functions
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
    𝒑₀ᵀ𝗠⁻¹, 𝒑₀ᵀ∂𝗠⁻¹∂x = cal∂𝗠∂x!(ap,x)
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

## set shape functions
function set𝝭(aps::Vector{T}) where T<:ReproducingKernel{SNode}
    for ap in aps
        set𝝭(ap)
    end
end
function set∇𝝭(aps::Vector{T}) where T<:ReproducingKernel{SNode}
    for ap in aps
        set∇𝝭(ap)
    end
end

function set𝝭(ap::ReproducingKernel{SNode})
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        i = ξ.id
        I = ξ.index[i]
        ξ̂ = Node(ξ)
        𝝭 = get𝝭(ap,ξ̂)
        for j in 1:length(𝓒)
            ξ.𝝭[:∂1][I+j] = 𝝭[j]
        end
    end
end

function set∇𝝭(ap::ReproducingKernel{SNode})
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        i = ξ.id
        I = ξ.index[i]
        ξ̂ = Node(ξ)
        𝝭,∂𝝭∂x,∂𝝭∂y,∂𝝭∂z = get∇𝝭(ap,ξ̂)
        for j in 1:length(𝓒)
            ξ.𝝭[:∂1][I+j] = 𝝭[j]
            ξ.𝝭[:∂x][I+j] = ∂𝝭∂x[j]
            ξ.𝝭[:∂y][I+j] = ∂𝝭∂y[j]
            ξ.𝝭[:∂z][I+j] = ∂𝝭∂z[j]
        end
    end
end

## convert
function Poi1(aps::Vector{T};renumbering::Bool=false) where T<:Approximator
    aps_ = Poi1[]
    𝓖 = Node[]
    if renumbering
        index, data = renumber(aps)
        for ap in aps
            i = ap.𝓒[1].id
            𝓒 = [Node(index[i],data)]
            push!(aps_,Poi1(𝓒,𝓖))
        end
    else
        for ap in aps
            𝓒 = [ap.𝓒[i]]
            push!(aps_,Poi1(𝓒,𝓖))
        end
    end
    return aps_
end
function Seg2(aps::Vector{T};renumbering::Bool=false) where T<:Approximator
    aps_ = Seg2[]
    𝓖 = Node[]
    if renumbering
        index, data = renumber(aps)
        for ap in aps
            i = ap.𝓒[1].id
            j = ap.𝓒[2].id
            𝓒 = [Node(index[i],data),Node(index[j],data)]
            push!(aps_,Seg2(𝓒,𝓖))
        end
        return aps_, data
    else
        for ap in aps
            𝓒 = [ap.𝓒[i] for i in 1:2]
            push!(aps_,Seg2(𝓒,𝓖))
        end
        return aps_
    end
end

function SegN{T,𝒑,𝑠,𝜙}(aps::Vector{A}) where {T<:AbstractNode,𝒑,𝑠,𝜙,A<:Approximator}
    aps_ = SegN{T,𝒑,𝑠,𝜙}[]
    𝗠 = Dict{Symbol,SymMat}()
    𝝭 = Dict{Symbol,Vector{Float64}}()
    𝓖 = T[]
    for ap in aps
        𝓒 = ap.𝓒
        push!(aps_,SegN{T,𝒑,𝑠,𝜙}(𝓒,𝗠,𝝭))
    end
    return aps_
end

## RK gradient smoothing
function set∇̃𝝭(gps::Vector{T},aps::Vector{S}) where{T<:ReproducingKernel,S<:ReproducingKernel}
    if length(gps) ≠ length(aps)
        error("Miss match element numbers")
    else
        for i in 1:length(gps)
            set∇̃𝝭(gps[i],aps[i])
        end
    end
end
function set∇̃𝝭(gp::SegN{SNode},ap::SegN{SNode})
    𝗚⁻¹ = cal𝗚!(gp)
    𝓒 = gp.𝓒
    𝓖 = gp.𝓖
    for ξ̂ in 𝓖
        𝒑̂ = get𝒑(gp,ξ̂)
        𝒑̂ᵀ𝗚⁻¹ = 𝒑̂*𝗚⁻¹
        ∂𝝭∂x = gp.𝝭[:∂x]
        fill!(∂𝝭∂x,0.0)
        for ξ in ap.𝓖
            w = ξ.w/2
            wᵇ = ξ.wᵇ
            n₁ = ξ.n₁
            𝝭 = get𝝭(ap,ξ)
            𝒑, ∂𝒑∂ξ = get∇𝒑(gp,ξ)
            W₁ = 𝒑̂ᵀ𝗚⁻¹*𝒑*n₁*wᵇ + 𝒑̂ᵀ𝗚⁻¹*∂𝒑∂ξ*n₁*w
            for i in 1:length(𝓒)
                ∂𝝭∂x[i] += 𝝭[i]*W₁
            end
        end
        for i in 1:length(𝓒)
            ξ̂.𝝭[:∂x][ξ̂.index[ξ̂.id]+i] = ∂𝝭∂x[i]
        end
    end
end

function set∇̃𝝭(gp::TriN{SNode},ap::TriN{SNode})
    𝗚⁻¹ = cal𝗚!(gp)
    𝓒 = gp.𝓒
    𝓖 = gp.𝓖
    for ξ̂ in 𝓖
        𝒑̂ = get𝒑(gp,ξ̂)
        𝒑̂ᵀ𝗚⁻¹ = 𝒑̂*𝗚⁻¹
        ∂𝝭∂x = gp.𝝭[:∂x]
        ∂𝝭∂y = gp.𝝭[:∂y]
        fill!(∂𝝭∂x,0.0)
        fill!(∂𝝭∂y,0.0)
        for ξ in ap.𝓖
            w = ξ.w
            wᵇ = ξ.wᵇ
            n₁ = ξ.n₁
            n₂ = ξ.n₂
            𝝭 = get𝝭(ap,ξ)
            𝒑, ∂𝒑∂ξ, ∂𝒑∂η = get∇𝒑(gp,ξ)
            b = 𝒑̂ᵀ𝗚⁻¹*∂𝒑∂ξ*n₁ + 𝒑̂ᵀ𝗚⁻¹*∂𝒑∂η*n₂
            W₁ = 𝒑̂ᵀ𝗚⁻¹*𝒑*n₁*wᵇ + b*w/2
            W₂ = 𝒑̂ᵀ𝗚⁻¹*𝒑*n₂*wᵇ + b*w/2
            for i in 1:length(𝓒)
                ∂𝝭∂x[i] += 𝝭[i]*W₁
                ∂𝝭∂y[i] += 𝝭[i]*W₂
            end
        end
        for i in 1:length(𝓒)
            ξ̂.𝝭[:∂x][ξ̂.index[ξ̂.id]+i] = ∂𝝭∂x[i]
            ξ̂.𝝭[:∂y][ξ̂.index[ξ̂.id]+i] = ∂𝝭∂y[i]
        end
    end
end

function set∇̃𝝭(gp::TetN{SNode},ap::TetN{SNode})
    𝗚⁻¹ = cal𝗚!(gp)
    𝓒 = gp.𝓒
    𝓖 = gp.𝓖
    for ξ̂ in 𝓖
        𝒑̂ = get𝒑(gp,ξ̂)
        𝒑̂ᵀ𝗚⁻¹ = 𝒑̂*𝗚⁻¹
        ∂𝝭∂x = gp.𝝭[:∂x]
        ∂𝝭∂y = gp.𝝭[:∂y]
        ∂𝝭∂z = gp.𝝭[:∂z]
        fill!(∂𝝭∂x,0.0)
        fill!(∂𝝭∂y,0.0)
        fill!(∂𝝭∂z,0.0)
        for ξ in ap.𝓖
            w = ξ.w
            wᵇ = ξ.wᵇ
            n₁ = ξ.n₁
            n₂ = ξ.n₂
            n₃ = ξ.n₃
            𝝭 = get𝝭(ap,ξ)
            𝒑, ∂𝒑∂ξ, ∂𝒑∂η, ∂𝒑∂γ = get∇𝒑(gp,ξ)
            b = 𝒑̂ᵀ𝗚⁻¹*∂𝒑∂ξ*n₁ + 𝒑̂ᵀ𝗚⁻¹*∂𝒑∂η*n₂ + 𝒑̂ᵀ𝗚⁻¹*∂𝒑∂γ*n₃
            W₁ = 𝒑̂ᵀ𝗚⁻¹*𝒑*n₁*wᵇ + b*w/3
            W₂ = 𝒑̂ᵀ𝗚⁻¹*𝒑*n₂*wᵇ + b*w/3
            W₃ = 𝒑̂ᵀ𝗚⁻¹*𝒑*n₃*wᵇ + b*w/3
            for i in 1:length(𝓒)
                ∂𝝭∂x[i] += 𝝭[i]*W₁
                ∂𝝭∂y[i] += 𝝭[i]*W₂
                ∂𝝭∂z[i] += 𝝭[i]*W₃
            end
        end
        for i in 1:length(𝓒)
            ξ̂.𝝭[:∂x][ξ̂.index[ξ̂.id]+i] = ∂𝝭∂x[i]
            ξ̂.𝝭[:∂y][ξ̂.index[ξ̂.id]+i] = ∂𝝭∂y[i]
            ξ̂.𝝭[:∂z][ξ̂.index[ξ̂.id]+i] = ∂𝝭∂z[i]
        end
    end
end

function renumber(aps::Vector{T}) where T<:Approximator
    index = Dict{Int,Int}()
    n = 0
    for ap in aps
        for x in ap.𝓒
            I = x.id
            if ~haskey(index,I)
                n += 1
                index[I] = n
            end
        end
    end
    data_ = aps[1].𝓒[1].data
    data = Dict(:x=>zeros(n),:y=>zeros(n),:z=>zeros(n))
    for (j,i) in index
        data[:x][i] = data_[:x][j]
        data[:y][i] = data_[:y][j]
        data[:z][i] = data_[:z][j]
    end
    return index, data
end
