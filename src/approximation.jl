
struct Element{T}<:AbstractElement{T}
    𝓒::Vector{Node}
    𝓖::Vector{Node}
end
Element{T}(𝓒::Vector{Node}) where T = Element{T}(𝓒,Node[])
function Element{T}(data::Dict{Symbol,Vector{Float64}},index::Int...) where T
    𝓒 = [Node(i,data) for i in index]
    𝓖 = Node[]
    return Element{T}(𝓒,𝓖)
end

## convert
Element{T}(a::S) where {T,S<:AbstractElement} = Element{T}(a.𝓒)
function Element{T}(as::Vector{S};renumbering::Bool=false) where {T,S<:AbstractElement}
    aps = Element{T}[]
    if renumbering
        index, data = renumber(aps)
        for a in as
            𝓒 = [Node(index[x.id],data) for x in a.𝓒]
            𝓖 = Node[]
            push!(aps,Element{T}(𝓒,𝓖))
        end
    else
        for a in as
            push!(aps,Element{T}(a))
        end
    end
    return aps
end

function Element{T}(a::AbstractElement,b::AbstractElement) where T
    𝓒 = a.𝓒
    𝓖 = get𝓖(a,b)
    𝓖 ≠ nothing ? Element{T}(𝓒,𝓖) : nothing
end
function Element{T}(as::Vector{A},bs::Vector{B}) where {T,A<:AbstractElement,B<:AbstractElement}
    aps = Element{T}[]
    for a in as
        for b in bs
            ap = Element{T}(a,b)
            ap ≠ nothing ? push!(aps,ap) : nothing
        end
    end
    return aps
end
## get𝒙
@inline get𝒙(ap::T,::Any) where T<:AbstractElement{:Poi1} = (ap.𝓒[1].x,ap.𝓒[1].y,ap.𝓒[1].z)
@inline get𝒙(ap::T,ξ::𝝃) where {T<:AbstractElement{:Seg2},𝝃<:AbstractNode} = get𝒙(ap,ξ.ξ)
@inline get𝒙(ap::T,ξ::𝝃) where {T<:AbstractElement{:Tri3},𝝃<:AbstractNode} = get𝒙(ap,ξ.ξ,ξ.η)

function get𝒙(ap::T,ξ::Float64) where T<:AbstractElement{:Seg2}
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
function get𝒙(ap::T,ξ::Float64,η::Float64) where T<:AbstractElement{:Tri3}
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

## get𝐽
@inline get𝐽(  ::T,::Any) where T<:AbstractElement{:Poi1} = 1.0
@inline get𝐽(ap::T,::Any) where T<:AbstractElement{:Seg2} = 0.5*get𝐿(ap)
@inline get𝐽(ap::T,::Any) where T<:AbstractElement{:Tri3} = get𝐴(ap)
## get𝑤
@inline get𝑤(  ::T,::Any) where T<:AbstractElement{:Poi1} = 1.0
@inline get𝑤(ap::T,ξ::𝝃) where {T<:AbstractElement{:Seg2},𝝃<:AbstractNode} = 0.5*get𝐿(ap)*ξ.w
@inline get𝑤(ap::T,ξ::𝝃) where {T<:AbstractElement{:Tri3},𝝃<:AbstractNode} = get𝐴(ap)*ξ.w
## get𝐿 get𝐴 get𝑉
@inline function get𝐿(ap::T) where T<:AbstractElement{:Seg2}
    x₁ = ap.𝓒[1].x
    y₁ = ap.𝓒[1].y
    z₁ = ap.𝓒[1].z
    x₂ = ap.𝓒[2].x
    y₂ = ap.𝓒[2].y
    z₂ = ap.𝓒[2].z
    return ((x₂-x₁)^2+(y₂-y₁)^2+(z₂-z₁)^2)^0.5
end
function get𝐴(ap::T) where T<:AbstractElement{:Tri3}
    x₁ = ap.𝓒[1].x
    y₁ = ap.𝓒[1].y
    z₁ = ap.𝓒[1].z
    x₂ = ap.𝓒[2].x
    y₂ = ap.𝓒[2].y
    z₂ = ap.𝓒[2].z
    x₃ = ap.𝓒[3].x
    y₃ = ap.𝓒[3].y
    z₃ = ap.𝓒[3].z
    𝐴₁ = 0.5*(y₁*z₂+y₂*z₃+y₃*z₁-y₂*z₁-y₃*z₂-y₁*z₃)
    𝐴₂ = 0.5*(z₁*x₂+z₂*x₃+z₃*x₁-z₂*x₁-z₃*x₂-z₁*x₃)
    𝐴₃ = 0.5*(x₁*y₂+x₂*y₃+x₃*y₁-x₂*y₁-x₃*y₂-x₁*y₃)
    return (𝐴₁^2 + 𝐴₂^2 + 𝐴₃^2)^0.5
end
## get𝒏
@inline get𝒏(ap::T) where T<:AbstractElement{:Seg2} = 1.0
@inline get𝒏(ap::T,ξ::𝝃) where {T<:AbstractElement{:Seg2},𝝃<:AbstractNode} = get𝒏(ap,ξ.ξ)
@inline get𝒏(ap::T,ξ::𝝃) where {T<:AbstractElement{:Tri3},𝝃<:AbstractNode} = get𝒏(ap,ξ.ξ,ξ.η)

function get𝒏(ap::T,ξ::Float64) where T<:AbstractElement{:Seg2}
    n₁ = 0.0
    n₁ += ξ == -1.0 ? -1.0 : 0.0
    n₁ += ξ ==  1.0 ?  1.0 : 0.0
    return n₁
end
function get𝒏(ap::T,ξ::Float64,η::Float64) where T<:AbstractElement{:Tri3}
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
## set𝒏!
function set𝒏!(aps::Vector{T}) where T<:AbstractElement
    for ap in aps
        set𝒏!(ap)
    end
end

function set𝒏!(ap::T) where T<:AbstractElement{:Seg2}
    𝓖 = ap.𝓖
    for ξ in 𝓖
        ξ.n₁ = get𝒏(ap,ξ)
    end
end

function set𝒏!(ap::T) where T<:AbstractElement{:Tri3}
    𝓖 = ap.𝓖
    for ξ in 𝓖
        ξ.n₁, ξ.n₂ = get𝒏(ap,ξ)
    end
end

## shape functions
# ------------- Poi1 ---------------
get𝝭(::Element{:Poi1},::Node) = 1.0
# ------------- Seg2 ---------------
@inline get𝝭(ap::Element{:Seg2},ξ::Node) = get𝝭(ap,ξ.ξ)
@inline get𝝭(ap::Element{:Seg2},ξ::Float64) = (0.5*(1-ξ),0.5*(1+ξ))
@inline function get∂𝝭∂x(ap::Element{:Seg2},::Any)
    𝐿 = get𝐿(ap)
    return (-1.0/𝐿,1.0/𝐿)
end
@inline get∂𝝭∂y(ap::Element{:Seg2},::Any) = (0.0,0.0)
@inline get∂𝝭∂z(ap::Element{:Seg2},::Any) = (0.0,0.0)
@inline get∇𝝭(ap::Element{:Seg2},ξ::Node) = (get𝝭(ap,ξ),get∂𝝭∂x(ap,ξ),(0.0,0.0),(0.0,0.0))
@inline function get∂𝝭∂𝑛(ap::Element{:Seg2},ξ::Node)
    n₁ = get𝒏(ap,ξ)
    𝐿 = get𝐿(ap)
    return (-n₁/𝐿,n₁/𝐿)
end
@inline get∇𝑛𝝭(ap::Element{:Seg2},ξ::Node) = (get𝝭(ap,ξ),get∂𝝭∂𝑛(ap,ξ))
# ------------- Tri3 ---------------

# ------------- Quad ---------------
get𝝭(ap::Element{:Quad},ξ::Node) = get𝝭(ap,ξ.ξ,ξ.η)
get∂𝝭∂ξ(ap::Element{:Quad},ξ::Node) = get∂𝝭∂ξ(ap,ξ.η)
get∂𝝭∂η(ap::Element{:Quad},ξ::Node) = get∂𝝭∂η(ap,ξ.ξ)

function get𝝭(ap::Element{:Quad},ξ::Float64,η::Float64)
    N₁ = 0.25*(1.0-ξ)*(1.0-η)
    N₂ = 0.25*(1.0+ξ)*(1.0-η)
    N₃ = 0.25*(1.0+ξ)*(1.0+η)
    N₄ = 0.25*(1.0-ξ)*(1.0+η)
    return (N₁,N₂,N₃,N₄)
end
function get∂𝝭∂ξ(ap::Element{:Quad},η::Float64)
    ∂N₁∂ξ = - 0.25*(1-η)
    ∂N₂∂ξ =   0.25*(1-η)
    ∂N₃∂ξ =   0.25*(1+η)
    ∂N₄∂ξ = - 0.25*(1+η)
    return (∂N₁∂ξ,∂N₂∂ξ,∂N₃∂ξ,∂N₄∂ξ)
end
function get∂𝝭∂η(ap::Element{:Quad},ξ::Float64)
    ∂N₁∂η = - 0.25*(1-ξ)
    ∂N₂∂η = - 0.25*(1+ξ)
    ∂N₃∂η =   0.25*(1+ξ)
    ∂N₄∂η =   0.25*(1-ξ)
    return (∂N₁∂η,∂N₂∂η,∂N₃∂η,∂N₄∂η)
end
function get∂𝝭∂x∂𝝭∂y(ap::Element{:Quad},ξ::Node)
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
get∇𝝭(ap::Element{:Quad},ξ::Node) = get𝝭(ap,ξ),get∂𝝭∂x∂𝝭∂y(ap,ξ)...,(0.0,0.0,0.0,0.0)
