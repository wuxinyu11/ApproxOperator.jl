"""
# Element
"""
struct Element{T}<:AbstractElement{T}
    𝓒::Vector{Node}
    𝓖::Vector{SNode}
end
Element{T}(𝓒::Vector{Node}) where T = Element{T}(𝓒,SNode[])
Element{T}(a::S) where {T,S<:AbstractElement} = Element{T}(a.𝓒)

"""
set𝒙!(ap::T,x::SNode) where T<:AbstractElement
get𝒙(ap::T,ξ::Float64...) where T<:AbstractElement
"""
function set𝒙!(aps::Vector{T}) where T<:AbstractElement
    nᵢ = getnᵢ(aps)
    data = getfield(aps[end].𝓖[end],:data)
    push!(data,:x=>(2,zeros(nᵢ)),:y=>(2,zeros(nᵢ)),:z=>(2,zeros(nᵢ)))
    set𝒙!.(aps)
end
function set𝒙!(ap::T) where T<:AbstractElement
    𝓖 = ap.𝓖
    for ξ in 𝓖
        x,y,z = get𝒙(ap,ξ)
        ξ.x = x
        ξ.y = y
        ξ.z = z
    end
end
@inline get𝒙(ap::T,::Any) where T<:AbstractElement{:Poi1} = (ap.𝓒[1].x,ap.𝓒[1].y,ap.𝓒[1].z)
@inline get𝒙(ap::T,ξ::SNode) where T<:AbstractElement{:Seg2} = get𝒙(ap,ξ.ξ)
@inline get𝒙(ap::T,ξ::SNode) where T<:AbstractElement{:Seg3} = get𝒙(ap,ξ.ξ)
@inline get𝒙(ap::T,ξ::SNode) where T<:AbstractElement{:Tri3} = get𝒙(ap,ξ.ξ,ξ.η)
@inline get𝒙(ap::T,ξ::SNode) where T<:AbstractElement{:Tri6} = get𝒙(ap,ξ.ξ,ξ.η)
@inline get𝒙(ap::T,ξ::SNode) where T<:AbstractElement{:Quad} = get𝒙(ap,ξ.ξ,ξ.η)

function get𝒙(ap::T,ξ::Float64) where T<:AbstractElement{:Seg2}
    x₁ = ap.𝓒[1].x
    y₁ = ap.𝓒[1].y
    z₁ = ap.𝓒[1].z
    x₂ = ap.𝓒[2].x
    y₂ = ap.𝓒[2].y
    z₂ = ap.𝓒[2].z
    N₁ = 0.5*(1.0-ξ)
    N₂ = 0.5*(1.0+ξ)
    return (x₁*N₁+x₂*N₂,y₁*N₁+y₂*N₂,z₁*N₁+z₂*N₂)
end
function get𝒙(ap::T,ξ::Float64) where T<:AbstractElement{:Seg3}
    x₁ = ap.𝓒[1].x
    y₁ = ap.𝓒[1].y
    z₁ = ap.𝓒[1].z
    x₂ = ap.𝓒[2].x
    y₂ = ap.𝓒[2].y
    z₂ = ap.𝓒[2].z
    x₃ = ap.𝓒[3].x
    y₃ = ap.𝓒[3].y
    z₃ = ap.𝓒[3].z
    N₁ = 0.5*ξ*(ξ-1.0)
    N₂ = 1.0-ξ^2
    N₃ = 0.5*ξ*(ξ+1.0)
    return (x₁*N₁+x₂*N₂+x₃*N₃,y₁*N₁+y₂*N₂+y₃*N₃,z₁*N₁+z₂*N₂+z₃*N₃)
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

function get𝒙(ap::T,ξ::Float64,η::Float64) where T<:AbstractElement{:Tri6}
    γ = 1.0-ξ-η
    x₁ = ap.𝓒[1].x;y₁ = ap.𝓒[1].y;z₁ = ap.𝓒[1].z
    x₂ = ap.𝓒[2].x;y₂ = ap.𝓒[2].y;z₂ = ap.𝓒[2].z
    x₃ = ap.𝓒[3].x;y₃ = ap.𝓒[3].y;z₃ = ap.𝓒[3].z
    x₄ = ap.𝓒[4].x;y₄ = ap.𝓒[4].y;z₄ = ap.𝓒[4].z
    x₅ = ap.𝓒[5].x;y₅ = ap.𝓒[5].y;z₅ = ap.𝓒[5].z
    x₆ = ap.𝓒[6].x;y₆ = ap.𝓒[6].y;z₆ = ap.𝓒[6].z
    N₁ = ξ*(2*ξ-1)
    N₂ = η*(2*η-1)
    N₃ = γ*(2*γ-1)
    N₄ = 4*ξ*η
    N₅ = 4*η*γ
    N₆ = 4*γ*ξ
    return (x₁*N₁+x₂*N₂+x₃*N₃+x₄*N₄+x₅*N₅+x₆*N₆,
            y₁*N₁+y₂*N₂+y₃*N₃+y₄*N₄+y₅*N₅+y₆*N₆,
            z₁*N₁+z₂*N₂+z₃*N₃+z₄*N₄+z₅*N₅+z₆*N₆)
end

function get𝒙(ap::T,ξ::Float64,η::Float64) where T<:AbstractElement{:Quad}
    x₁ = ap.𝓒[1].x
    y₁ = ap.𝓒[1].y
    z₁ = ap.𝓒[1].z
    x₂ = ap.𝓒[2].x
    y₂ = ap.𝓒[2].y
    z₂ = ap.𝓒[2].z
    x₃ = ap.𝓒[3].x
    y₃ = ap.𝓒[3].y
    z₃ = ap.𝓒[3].z
    x₄ = ap.𝓒[4].x
    y₄ = ap.𝓒[4].y
    z₄ = ap.𝓒[4].z
    N₁,N₂,N₃,N₄ = get𝝭(ap,ξ,η)
    return (x₁*N₁+x₂*N₂+x₃*N₃+x₄*N₄,y₁*N₁+y₂*N₂+y₃*N₃+y₄*N₄,z₁*N₁+z₂*N₂+z₃*N₃+z₄*N₄)
end

function get𝑱(ap::T,ξ::SNode) where T<:AbstractElement{:Quad}
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

"""
get𝐽(ap::T,x::SNode) where T<:AbstractElement
"""
@inline get𝐽(  ::T,::Any) where T<:AbstractElement{:Poi1} = 1.0
@inline get𝐽(ap::T,::Any) where T<:AbstractElement{:Seg2} = 0.5*get𝐿(ap)
@inline get𝐽(ap::T,::Any) where T<:AbstractElement{:Tri3} = 2.0*get𝐴(ap)
@inline function get𝐽(ap::T,ξ::SNode) where T<:AbstractElement{:Quad}
    J₁₁,J₂₁,J₁₂,J₂₂ = get𝑱(ap,ξ)
    return J₁₁*J₂₂-J₂₁*J₁₂
end

"""
set𝑤!(ap::T) where T<:AbstractElement
get𝑤(ap::T,x::SNode) where T<:AbstractElement
"""
function set𝑤!(aps::Vector{T}) where T<:AbstractElement
    nᵢ = getnᵢ(aps)
    data = getfield(aps[end].𝓖[end],:data)
    push!(data,:𝑤=>(2,zeros(nᵢ)))
    set𝑤!.(aps)
end
function set𝑤!(ap::T) where T<:AbstractElement
    𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = get𝑤(ap,ξ)
        ξ.𝑤 = 𝑤
    end
end
@inline get𝑤(  ::T,::Any) where T<:AbstractElement{:Poi1} = 1.0
@inline get𝑤(ap::T,ξ::SNode) where T<:AbstractElement{:Seg2} = 0.5*get𝐿(ap)*ξ.w
@inline get𝑤(ap::T,ξ::SNode) where T<:AbstractElement{:Seg3} = 0.5*get𝐿(ap)*ξ.w
@inline get𝑤(ap::T,ξ::SNode) where T<:AbstractElement{:Tri3} = get𝐴(ap)*ξ.w
@inline get𝑤(ap::T,ξ::SNode) where T<:AbstractElement{:Tri6} = get𝐴(ap)*ξ.w
@inline get𝑤(ap::T,ξ::SNode) where T<:AbstractElement{:Quad} = get𝐽(ap,ξ)*ξ.w

"""
get𝐿,get𝐴,get𝑉
"""
function set𝐿!(aps::Vector{T}) where T<:AbstractElement
    nₑ = length(aps)
    push!(getfield(aps[1].𝓖[1],:data),:𝐿=>(3,zeros(nₑ)))
    set𝐿!.(aps)
end
function set𝐿!(ap::T) where T<:AbstractElement
    𝐿 = get𝐿(ap)
    ap.𝓖[1].𝐿 = 𝐿
end
function set𝐴!(aps::Vector{T}) where T<:AbstractElement
    nₑ = length(aps)
    push!(getfield(aps[1].𝓖[1],:data),:𝐴=>(3,zeros(nₑ)))
    set𝐴!.(aps)
end
function set𝐴!(ap::T) where T<:AbstractElement
    𝐴 = get𝐴(ap)
    ap.𝓖[1].𝐴 = 𝐴
end
function set𝑉!(aps::Vector{T}) where T<:AbstractElement
    nₑ = length(aps)
    push!(getfield(aps[1].𝓖[1],:data),:𝑉=>(3,zeros(nₑ)))
    set𝑉!.(aps)
end
function set𝒙ₘ!(aps::Vector{T}) where T<:AbstractElement
    nₑ = length(aps)
    push!(getfield(aps[1].𝓖[1],:data),:xₘ=>(3,zeros(nₑ)),:yₘ=>(3,zeros(nₑ)))
    set𝒙ₘ!.(aps)
end
function set𝒙ₘ!(ap::T) where T<:AbstractElement
    xₘ,yₘ = get𝒙ₘ(ap)
    ap.𝓖[1].xₘ = xₘ
    ap.𝓖[1].yₘ = yₘ
end
function setm2!(aps::Vector{T}) where T<:AbstractElement
    nₑ = length(aps)
    push!(getfield(aps[1].𝓖[1],:data),:m₂₀=>(3,zeros(nₑ)),:m₁₁=>(3,zeros(nₑ)),:m₀₂=>(3,zeros(nₑ)))
    setm2!.(aps)
end
function setm2!(ap::T) where T<:AbstractElement
    m₂₀,m₁₁,m₀₂ = getm2(ap)
    ap.𝓖[1].m₂₀ = m₂₀
    ap.𝓖[1].m₁₁ = m₁₁
    ap.𝓖[1].m₀₂ = m₀₂
end

function get𝐿(ap::T) where T<:AbstractElement{:Seg2}
    x₁ = ap.𝓒[1].x
    y₁ = ap.𝓒[1].y
    z₁ = ap.𝓒[1].z
    x₂ = ap.𝓒[2].x
    y₂ = ap.𝓒[2].y
    z₂ = ap.𝓒[2].z
    return ((x₂-x₁)^2+(y₂-y₁)^2+(z₂-z₁)^2)^0.5
end
function get𝐿(ap::T) where T<:AbstractElement{:Seg3}
    x₁ = ap.𝓒[1].x
    y₁ = ap.𝓒[1].y
    z₁ = ap.𝓒[1].z
    x₂ = ap.𝓒[3].x
    y₂ = ap.𝓒[3].y
    z₂ = ap.𝓒[3].z
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
function get𝐴(ap::T) where T<:AbstractElement{:Tri6}
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
function get𝐴(ap::T) where T<:AbstractElement{:Vor2}
    𝓒 = ap.𝓒
    nᵥ = length(𝓒)
    𝐴 = 0.0
    for i in 1:nᵥ-1
        x₁ = 𝓒[i].x
        y₁ = 𝓒[i].y
        x₂ = 𝓒[i+1].x
        y₂ = 𝓒[i+1].y
        𝐴 += x₁*y₂-x₂*y₁
    end
    x₁ = 𝓒[nᵥ].x
    y₁ = 𝓒[nᵥ].y
    x₂ = 𝓒[1].x
    y₂ = 𝓒[1].y
    𝐴 += x₁*y₂-x₂*y₁
    𝐴 *= 0.5
    return 𝐴
end
function get𝐴(ap::T) where T<:AbstractElement{:Vor2}
    𝓒 = ap.𝓒
    nᵥ = length(𝓒)
    𝐴 = 0.0
    for i in 1:nᵥ
        𝐴 += i≠nᵥ ? 𝓒[i].x*𝓒[i+1].y-𝓒[i+1].x*𝓒[i].y : 𝓒[i].x*𝓒[1].y-𝓒[1].x*𝓒[i].y
    end
    return 0.5*𝐴
end

function get𝒙ₘ(ap::AbstractElement)
    𝓒 = ap.𝓒
    nᵥ = length(𝓒)
    𝐴 = ap.𝓖[1].𝐴
    xₘ = 0.0
    yₘ = 0.0
    for i in 1:nᵥ-1
        x₁ = 𝓒[i].x
        y₁ = 𝓒[i].y
        x₂ = 𝓒[i+1].x
        y₂ = 𝓒[i+1].y
        xₘ += (x₁*y₂-x₂*y₁)*(x₁+x₂)
        yₘ += (x₁*y₂-x₂*y₁)*(y₁+y₂)
    end
    x₁ = 𝓒[nᵥ].x
    y₁ = 𝓒[nᵥ].y
    x₂ = 𝓒[1].x
    y₂ = 𝓒[1].y
    xₘ += (x₁*y₂-x₂*y₁)*(x₁+x₂)
    yₘ += (x₁*y₂-x₂*y₁)*(y₁+y₂)
    xₘ /= 6.0*𝐴
    yₘ /= 6.0*𝐴
    return xₘ,yₘ
end
function getm2(ap::AbstractElement)
    𝓒 = ap.𝓒
    nᵥ = length(𝓒)
    𝐴 = ap.𝓖[1].𝐴
    xₘ = ap.𝓖[1].xₘ
    yₘ = ap.𝓖[1].yₘ
    m₂₀ = 0.0
    m₁₁ = 0.0
    m₀₂ = 0.0
    for i in 1:nᵥ-1
        x₁ = 𝓒[i].x
        y₁ = 𝓒[i].y
        x₂ = 𝓒[i+1].x
        y₂ = 𝓒[i+1].y
        m₂₀ += (x₁*y₂-x₂*y₁)*(x₁^2+x₁*x₂+x₂^2)
        m₁₁ += (x₁*y₂-x₂*y₁)*(2*x₁*y₁+x₁*y₂+x₂*y₁+2*x₂*y₂)
        m₀₂ += (x₁*y₂-x₂*y₁)*(y₁^2+y₁*y₂+y₂^2)
    end
    x₁ = 𝓒[nᵥ].x
    y₁ = 𝓒[nᵥ].y
    x₂ = 𝓒[1].x
    y₂ = 𝓒[1].y
    m₂₀ += (x₁*y₂-x₂*y₁)*(x₁^2+x₁*x₂+x₂^2)
    m₁₁ += (x₁*y₂-x₂*y₁)*(2*x₁*y₁+x₁*y₂+x₂*y₁+2*x₂*y₂)
    m₀₂ += (x₁*y₂-x₂*y₁)*(y₁^2+y₁*y₂+y₂^2)
    m₂₀ /= 12.0*𝐴
    m₁₁ /= 24.0*𝐴
    m₀₂ /= 12.0*𝐴
    return m₂₀,m₁₁,m₀₂
end

"""
set𝒏!
"""
function set𝒏!(aps::Vector{T}) where T<:AbstractElement{:Seg2}
    data = getfield(aps[1].𝓖[1],:data)
    n = length(aps)
    push!(data,:n₁=>(3,zeros(n)))
    push!(data,:n₂=>(3,zeros(n)))
    for ap in aps
        set𝒏!(ap)
    end
end

function set𝒏!(ap::T) where T<:AbstractElement{:Seg2}
    x₁ = ap.𝓒[1].x
    y₁ = ap.𝓒[1].y
    x₂ = ap.𝓒[2].x
    y₂ = ap.𝓒[2].y
    𝐿 = get𝐿(ap)
    n₁ = (y₂-y₁)/𝐿
    n₂ = (x₁-x₂)/𝐿
    ap.𝓖[1].n₁ = n₁
    ap.𝓖[1].n₂ = n₂
end

"""
set𝐷!
get𝐷!
"""
function set𝑫!(aps::Vector{T}) where T<:AbstractElement{:Seg2}
    n = getnᵢ(aps)
    data = getfield(aps[1].𝓖[1],:data)
    push!(data,:D₁=>(2,zeros(n)))
    for ap in aps
        set𝑫!(ap)
    end
end

function set𝑫!(ap::T) where T<:AbstractElement{:Seg2}
    for ξ in ap.𝓖
        ξ.D₁ = ξ.ξ == -1.0 ? -1.0 : 0.0
        ξ.D₁ = ξ.ξ ==  1.0 ?  1.0 : 0.0
    end
end

function set𝑫!(aps::Vector{T}) where T<:AbstractElement{:Tri3}
    data = getfield(aps[1].𝓖[1],:data)
    n = getnᵢ(aps)
    nₑ = length(aps)
    push!(data,:D₁=>(2,zeros(n)))
    push!(data,:D₂=>(2,zeros(n)))
    push!(data,:D₁₁=>(3,zeros(nₑ)))
    push!(data,:D₁₂=>(3,zeros(nₑ)))
    push!(data,:D₂₁=>(3,zeros(nₑ)))
    push!(data,:D₂₂=>(3,zeros(nₑ)))
    push!(data,:D₃₁=>(3,zeros(nₑ)))
    push!(data,:D₃₂=>(3,zeros(nₑ)))
    for ap in aps
        set𝑫!(ap)
    end
end

function set𝑫!(ap::T) where T<:AbstractElement{:Tri3}
    x₁ = ap.𝓒[1].x
    y₁ = ap.𝓒[1].y
    x₂ = ap.𝓒[2].x
    y₂ = ap.𝓒[2].y
    x₃ = ap.𝓒[3].x
    y₃ = ap.𝓒[3].y
    D₁₁ = y₃-y₂
    D₂₁ = y₁-y₃
    D₃₁ = y₂-y₁
    D₁₂ = x₂-x₃
    D₂₂ = x₃-x₁
    D₃₂ = x₁-x₂
    ap.𝓖[1].D₁₁ = D₁₁
    ap.𝓖[1].D₂₁ = D₂₁
    ap.𝓖[1].D₃₁ = D₃₁
    ap.𝓖[1].D₁₂ = D₁₂
    ap.𝓖[1].D₂₂ = D₂₂
    ap.𝓖[1].D₃₂ = D₃₂
    for ξ in ap.𝓖
        if ξ.ξ ≈ 0.0 (ξ.D₁ += D₁₁;ξ.D₂ += D₁₂) end
        if ξ.η ≈ 0.0 (ξ.D₁ += D₂₁;ξ.D₂ += D₂₂) end 
        if ξ.ξ+ξ.η ≈ 1.0 (ξ.D₁ += D₃₁;ξ.D₂ += D₃₂) end
    end
end

"""
setgeometry!(ap::T) where T<:AbstractElement
"""
function setgeometry!(aps::Vector{T}) where T<:AbstractElement
    set𝒙!(aps)
    set𝑤!(aps)
    if T<:AbstractElement{:Seg2}
        set𝐿!(aps)
    elseif T<:AbstractElement{:Tri3}
        set𝐴!(aps)
    elseif T<:AbstractElement{:Tet4}
        set𝑉!(aps)
    end
end

"""
set𝝭!
"""
function set𝝭!(aps::Vector{T}) where T<:AbstractElement
    for ap in aps
        set𝝭!(ap)
    end
end

function set𝝭!(ap::Element{S}) where S
    𝓖 = ap.𝓖
    for ξ in 𝓖
        N = get𝝭(ap,ξ)
        for i in 1:length(ap.𝓒)
            𝝭 = ξ[:𝝭]
            𝝭[i] = N[i]
        end
    end
end

"""
get𝝭(ap::Element,ξ::SNode)
"""
# ------------- Poi1 ---------------
function set𝝭!(ap::Element{:Poi1},x::SNode)
    𝝭 = x[:𝝭]
    𝝭[1] = 1.0
end

# ------------- Seg2 ---------------
function set𝝭!(ap::Element{:Seg2},x::SNode)
    𝝭 = x[:𝝭]
    𝝭[1] = 0.5*(1.0-x.ξ)
    𝝭[2] = 0.5*(1.0+x.ξ)
end

function set∇𝝭!(ap::Element{:Seg2},x::SNode)
    𝐿 = get𝐿(ap)
    ∂𝝭∂x = x[:∂𝝭∂x]
    ∂𝝭∂x[1] = -1.0/𝐿
    ∂𝝭∂x[2] = 1.0/𝐿
end

# ------------- Seg3 ---------------
function set𝝭!(ap::Element{:Seg3},x::SNode)
    𝝭 = x[:𝝭]
    ξ = x.ξ
    𝝭[1] = 0.5*ξ*(ξ-1.0)
    𝝭[2] = 1.0-ξ^2
    𝝭[3] = 0.5*ξ*(ξ+1.0)
end

function set∇𝝭!(ap::Element{:Seg3},x::SNode)
    𝐿 = get𝐿(ap)
    ∂𝝭∂x = x[:∂𝝭∂x]
    x₁ = ap.𝓒[1].x
    x₂ = ap.𝓒[2].x
    x₃ = ap.𝓒[3].x
    ξ = x.ξ
    ∂𝝭∂x[1] = (ξ-0.5)*2/𝐿
    ∂𝝭∂x[2] = -2.0*ξ*2/𝐿
    ∂𝝭∂x[3] = (ξ+0.5)*2/𝐿
end


# ------------- Tri3 ---------------
function set𝝭!(ap::Element{:Tri3},x::SNode)
    𝝭 = x[:𝝭]
    𝝭[1] = x.ξ
    𝝭[2] = x.η
    𝝭[3] = 1.0-x.ξ-x.η
end
function set∇𝝭!(ap::Element{:Tri3},x::SNode)
    𝐴 = get𝐴(ap)
    x₁ = ap.𝓒[1].x
    x₂ = ap.𝓒[2].x
    x₃ = ap.𝓒[3].x
    y₁ = ap.𝓒[1].y
    y₂ = ap.𝓒[2].y
    y₃ = ap.𝓒[3].y
    ∂𝝭∂x = x[:∂𝝭∂x]
    ∂𝝭∂y = x[:∂𝝭∂y]
    ∂𝝭∂x[1] = (y₂-y₃)/2.0/𝐴
    ∂𝝭∂x[2] = (y₃-y₁)/2.0/𝐴
    ∂𝝭∂x[3] = (y₁-y₂)/2.0/𝐴
    ∂𝝭∂y[1] = (x₃-x₂)/2.0/𝐴
    ∂𝝭∂y[2] = (x₁-x₃)/2.0/𝐴
    ∂𝝭∂y[3] = (x₂-x₁)/2.0/𝐴
end

# ------------- Tri6 ---------------
function set𝝭!(ap::Element{:Tri6},x::SNode)
    𝝭 = x[:𝝭]
    ξ = x.ξ
    η = x.η
    γ = 1.0-ξ-η
    𝝭[1] = ξ*(2*ξ-1)
    𝝭[2] = η*(2*η-1)
    𝝭[3] = γ*(2*γ-1)
    𝝭[4] = 4*ξ*η
    𝝭[5] = 4*η*γ
    𝝭[6] = 4*γ*ξ
end
function set∇𝝭!(ap::Element{:Tri6},x::SNode)
    𝐴 = get𝐴(ap)
    ξ = x.ξ
    η = x.η
    γ = 1.0-ξ-η
    x₁ = ap.𝓒[1].x
    x₂ = ap.𝓒[2].x
    x₃ = ap.𝓒[3].x
    y₁ = ap.𝓒[1].y
    y₂ = ap.𝓒[2].y
    y₃ = ap.𝓒[3].y
    J₁₁ = (y₂-y₃)/2.0/𝐴
    J₂₁ = (y₃-y₁)/2.0/𝐴
    J₃₁ = (y₁-y₂)/2.0/𝐴
    J₁₂ = (x₃-x₂)/2.0/𝐴
    J₂₂ = (x₁-x₃)/2.0/𝐴
    J₃₂ = (x₂-x₁)/2.0/𝐴
    ∂𝝭∂x = x[:∂𝝭∂x]
    ∂𝝭∂y = x[:∂𝝭∂y]
    ∂𝝭∂x[1] = (4.0*ξ-1.0)*J₁₁
    ∂𝝭∂x[2] = (4.0*η-1.0)*J₂₁
    ∂𝝭∂x[3] = (4.0*γ-1.0)*J₃₁
    ∂𝝭∂x[4] = 4.0*(η*J₁₁+ξ*J₂₁)
    ∂𝝭∂x[5] = 4.0*(γ*J₂₁+η*J₃₁)
    ∂𝝭∂x[6] = 4.0*(ξ*J₃₁+γ*J₁₁)
    ∂𝝭∂y[1] = (4.0*ξ-1.0)*J₁₂
    ∂𝝭∂y[2] = (4.0*η-1.0)*J₂₂
    ∂𝝭∂y[3] = (4.0*γ-1.0)*J₃₂
    ∂𝝭∂y[4] = 4.0*(η*J₁₂+ξ*J₂₂)
    ∂𝝭∂y[5] = 4.0*(γ*J₂₂+η*J₃₂)
    ∂𝝭∂y[6] = 4.0*(ξ*J₃₂+γ*J₁₂)
end

# ------------- Quad ---------------
function set𝝭!(ap::Element{:Quad},x::SNode)
    ξ = x.ξ
    η = x.η
    𝝭 = x[:𝝭]
    𝝭[1] = 0.25*(1.0-ξ)*(1.0-η)
    𝝭[2] = 0.25*(1.0+ξ)*(1.0-η)
    𝝭[3] = 0.25*(1.0+ξ)*(1.0+η)
    𝝭[4] = 0.25*(1.0-ξ)*(1.0+η)
end

function set∇𝝭!(ap::Element{:Quad},x::SNode)
    x₁ = ap.𝓒[1].x
    x₂ = ap.𝓒[2].x
    x₃ = ap.𝓒[3].x
    x₄ = ap.𝓒[4].x
    y₁ = ap.𝓒[1].y
    y₂ = ap.𝓒[2].y
    y₃ = ap.𝓒[3].y
    y₄ = ap.𝓒[4].y
    ∂N₁∂ξ,∂N₂∂ξ,∂N₃∂ξ,∂N₄∂ξ = get∂𝝭∂ξ(ap,x)
    ∂N₁∂η,∂N₂∂η,∂N₃∂η,∂N₄∂η = get∂𝝭∂η(ap,x)
    ∂x∂ξ = ∂N₁∂ξ*x₁ + ∂N₂∂ξ*x₂ + ∂N₃∂ξ*x₃ + ∂N₄∂ξ*x₄
    ∂x∂η = ∂N₁∂η*x₁ + ∂N₂∂η*x₂ + ∂N₃∂η*x₃ + ∂N₄∂η*x₄
    ∂y∂ξ = ∂N₁∂ξ*y₁ + ∂N₂∂ξ*y₂ + ∂N₃∂ξ*y₃ + ∂N₄∂ξ*y₄
    ∂y∂η = ∂N₁∂η*y₁ + ∂N₂∂η*y₂ + ∂N₃∂η*y₃ + ∂N₄∂η*y₄
    detJ = ∂x∂ξ*∂y∂η - ∂x∂η*∂y∂ξ
    ∂ξ∂x =   ∂y∂η/detJ
    ∂η∂x = - ∂y∂ξ/detJ
    ∂ξ∂y = - ∂x∂η/detJ
    ∂η∂y =   ∂x∂ξ/detJ
    ∂𝝭∂x = x[:∂𝝭∂x]
    ∂𝝭∂y = x[:∂𝝭∂y]
    ∂𝝭∂x[1] = ∂N₁∂ξ*∂ξ∂x + ∂N₁∂η*∂η∂x
    ∂𝝭∂x[2] = ∂N₂∂ξ*∂ξ∂x + ∂N₂∂η*∂η∂x
    ∂𝝭∂x[3] = ∂N₃∂ξ*∂ξ∂x + ∂N₃∂η*∂η∂x
    ∂𝝭∂x[4] = ∂N₄∂ξ*∂ξ∂x + ∂N₄∂η*∂η∂x
    ∂𝝭∂y[1] = ∂N₁∂ξ*∂ξ∂y + ∂N₁∂η*∂η∂y
    ∂𝝭∂y[2] = ∂N₂∂ξ*∂ξ∂y + ∂N₂∂η*∂η∂y
    ∂𝝭∂y[3] = ∂N₃∂ξ*∂ξ∂y + ∂N₃∂η*∂η∂y
    ∂𝝭∂y[4] = ∂N₄∂ξ*∂ξ∂y + ∂N₄∂η*∂η∂y
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

"""
Crouzeix-Raviart element
"""
struct TRElement{T}<:AbstractElement{T}
    𝓒::Vector{GNode}
    𝓖::Vector{SNode}
end

function set𝝭!(ap::TRElement{:Tri3},x::SNode)
    ξ₁ = x.ξ
    ξ₂ = x.η
    ξ₃ = 1.0-x.ξ-x.η
    N₁ = ξ₂+ξ₃-ξ₁
    N₂ = ξ₃+ξ₁-ξ₂
    N₃ = ξ₁+ξ₂-ξ₃
    𝝭 = x[:𝝭]
    𝝭[1] = N₁
    𝝭[2] = N₂
    𝝭[3] = N₃
end

function set∇𝝭!(ap::TRElement{:Tri3},x::SNode)
    x₁ = ap.𝓒[1].x
    x₂ = ap.𝓒[2].x
    x₃ = ap.𝓒[3].x
    y₁ = ap.𝓒[1].y
    y₂ = ap.𝓒[2].y
    y₃ = ap.𝓒[3].y
    𝐴 = get𝐴(ap)
    ∂𝝭∂x = x[:∂𝝭∂x]
    ∂𝝭∂y = x[:∂𝝭∂y]
    ∂𝝭∂x[1] = (y₃-y₂)/𝐴
    ∂𝝭∂x[2] = (y₁-y₃)/𝐴
    ∂𝝭∂x[3] = (y₂-y₁)/𝐴
    ∂𝝭∂y[1] = (x₂-x₃)/𝐴
    ∂𝝭∂y[2] = (x₃-x₁)/𝐴
    ∂𝝭∂y[3] = (x₁-x₂)/𝐴
end

function set∇̃𝝭!(ap::TRElement{:Tri3},x::SNode)
    x₁ = ap.𝓒[1].x
    x₂ = ap.𝓒[2].x
    x₃ = ap.𝓒[3].x
    y₁ = ap.𝓒[1].y
    y₂ = ap.𝓒[2].y
    y₃ = ap.𝓒[3].y
    𝐴 = get𝐴(ap)
    ∂𝝭∂x = x[:∂𝝭∂x]
    ∂𝝭∂y = x[:∂𝝭∂y]
    D₁₁ =  y₃-y₂
    D₁₂ =  x₂-x₃
    D₂₁ =  y₁-y₃
    D₂₂ =  x₃-x₁
    D₃₁ =  y₂-y₁
    D₃₂ =  x₁-x₂
    N = x[:𝝭]
    ∂𝝭∂x[1] = (D₁₁*N[1]-D₂₁*N[2]-D₃₁*N[3])/𝐴
    ∂𝝭∂x[2] = (D₂₁*N[1]-D₃₁*N[2]-D₁₁*N[3])/𝐴
    ∂𝝭∂x[3] = (D₃₁*N[1]-D₁₁*N[2]-D₂₁*N[3])/𝐴
    ∂𝝭∂y[1] = (D₁₂*N[1]-D₂₂*N[2]-D₃₂*N[3])/𝐴
    ∂𝝭∂y[2] = (D₂₂*N[1]-D₃₂*N[2]-D₁₂*N[3])/𝐴
    ∂𝝭∂y[3] = (D₃₂*N[1]-D₁₂*N[2]-D₂₂*N[3])/𝐴
end

for set𝝭 in (:set𝝭!,:set∇𝝭!,:set∇̄𝝭!,:set𝝭̄!,:set∇̃𝝭!)
    @eval begin
        function $set𝝭(aps::Vector{T}) where T<:AbstractElement
            for ap in aps
                𝓖 = ap.𝓖
                for 𝒙 in 𝓖
                    $set𝝭(ap,𝒙)
                end
            end
        end
    end
end
"""
⊆,∩
"""
function issubset(a::T,b::S) where {T<:AbstractElement{:Poi1},S<:AbstractElement{:Seg2}}
    i = findfirst(x->x==a.𝓒[1],b.𝓒)
    return i ≠ nothing && i ≤ 2
end

@inline intersect(a::T,b::T) where T<:AbstractElement = a.𝓒 == b.𝓒 ? a : nothing
@inline function intersect(a::T,b::S) where {T<:AbstractElement{:Seg2},S<:AbstractElement{:Poi1}}
    i = findfirst(x->x==b.𝓒[1],a.𝓒)
    return i ≠ nothing && i ≤ 2 ? a : nothing
end
@inline function intersect(a::T,b::S) where {T<:AbstractElement{:Tri3},S<:AbstractElement{:Poi1}}
    i = findfirst(x->x==b.𝓒[1],a.𝓒)
    return i ≠ nothing && i ≤ 3 ? a : nothing
end
@inline function intersect(a::T,b::S) where {T<:AbstractElement,S<:AbstractElement{:Seg2}}
    i = findfirst(x->x==b.𝓒[1],a.𝓒)
    j = findfirst(x->x==b.𝓒[2],a.𝓒)
    return i ≠ nothing && j ≠ nothing ? a : nothing
end
@inline function intersect(a::T,b::S) where {T<:AbstractElement{:Tri3},S<:AbstractElement{:Seg2}}
    i = findfirst(x->x==b.𝓒[1],a.𝓒)
    j = findfirst(x->x==b.𝓒[2],a.𝓒)
    return i ≠ nothing && j ≠ nothing && i ≤ 3 && j ≤ 3 ? a : nothing
end
@inline function intersect(a::T,b::S) where {T<:TRElement{:Tri3},S<:AbstractElement{:Seg2}}
    i = findfirst(x->x.𝑖==b.𝓒[1].𝐼, a.𝓒)
    j = findfirst(x->x.𝑖==b.𝓒[2].𝐼, a.𝓒)
    return i ≠ nothing && j ≠ nothing && i ≤ 3 && j ≤ 3 ? a : nothing
end
function intersect(as::Vector{T},bs::Vector{S}) where {T<:AbstractElement,S<:AbstractElement}
    aps = T[]
    for b in bs
        for a in as
            ap = a∩b
            ap ≠ nothing ? push!(aps,ap) : nothing
        end
    end
    return aps
end

function intersect(as::Vector{T},bs::Vector{S}) where {T<:AbstractElement{:Vor2},S<:AbstractElement{:Seg2}}
    aps = T[]
    ids = Int[]
    for b in bs
        for (i,a) in enumerate(as)
            ap = a∩b
            ap ≠ nothing ? (push!(aps,ap);push!(ids,i)) : nothing
        end
    end
    return aps, ids
end

"""
getnₚ,getnᵢ,getnₛ
"""
getnₚ(ap::T) where T<:AbstractElement = length(getfield(ap.𝓒[1],:data)[:x][2])
@inline getnₚ(aps::Vector{T}) where T<:AbstractElement = getnₚ(aps[1])
function getnₚ(aps::Vector{T}) where T<:TRElement
    nₚ = 0
    for ap in aps
        for x in ap.𝓒
            nₚ = max(nₚ,x.𝐼)
        end
    end
    return nₚ
end

function getnᵢ(aps::Vector{T}) where T<:AbstractElement
    ap = aps[end]
    ξ = ap.𝓖[end]
    return ξ.𝐺
end

function getnₛ(aps::Vector{T}) where T<:AbstractElement
    ap = aps[end]
    ξ = ap.𝓖[end]
    return ξ.𝑠 + length(ap.𝓒)
end

