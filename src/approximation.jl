
##
@inline function getx(ap::A,ξ::N) where {A<:Approximator,N<:AbstractNode}
    x = 0.0
    y = 0.0
    z = 0.0
    N = get𝝭(ap,ξ)
    for i in 1:length(𝓒)
        x += ap.𝓒[i].x*N[i]
        y += ap.𝓒[i].y*N[i]
        z += ap.𝓒[i].z*N[i]
    end
    return (x,y,z)
end

## AbstractSeg
@inline getw(ap::Seg2,ξ::T) where T<:AbstractNode = 0.5*ap.L*ξ.w
# ---------------- Seg2 -------------------
struct Seg2{T}<:Approximator where T<:AbstractNode
    𝓒::Vector{Node}
    𝓖::Vector{T}
    𝓡::Vector{SamplePoint}
    L::Float64
end

# constructions of Seg2
function Seg2(𝓒::Vector{Node},𝓖::Vector{T}) where T<:AbstractNode
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
