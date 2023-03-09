
@inline get𝒙(ap::T,ξ::Node) where T<:AbstractElement{:Seg3} = get𝒙(ap,ξ.ξ)

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

@inline get𝑤(ap::T,ξ::Node) where T<:AbstractElement{:Seg3} = 0.5*get𝐿(ap)*ξ.w

function get𝐿(ap::T) where T<:AbstractElement{:Seg3}
    x₁ = ap.𝓒[1].x
    y₁ = ap.𝓒[1].y
    z₁ = ap.𝓒[1].z
    x₂ = ap.𝓒[3].x
    y₂ = ap.𝓒[3].y
    z₂ = ap.𝓒[3].z
    return ((x₂-x₁)^2+(y₂-y₁)^2+(z₂-z₁)^2)^0.5
end

function set𝝭!(ap::Element{:Seg3},x::Node)
    𝝭 = x[:𝝭]
    ξ = x.ξ
    𝝭[1] = 0.5*ξ*(ξ-1.0)
    𝝭[2] = 1.0-ξ^2
    𝝭[3] = 0.5*ξ*(ξ+1.0)
end

function set∇𝝭!(ap::Element{:Seg3},x::Node)
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