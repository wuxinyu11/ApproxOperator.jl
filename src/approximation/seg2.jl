
@inline get𝒙(ap::T,ξ::Node) where T<:AbstractElement{:Seg2} = get𝒙(ap,ξ.ξ)

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

@inline get𝐽(ap::T,::Any) where T<:AbstractElement{:Seg2} = 0.5*get𝐿(ap)

@inline get𝑤(ap::T,ξ::Node) where T<:AbstractElement{:Seg2} = 0.5*get𝐿(ap)*ξ.w

function get𝐿(ap::T) where T<:AbstractElement{:Seg2}
    x₁ = ap.𝓒[1].x
    y₁ = ap.𝓒[1].y
    z₁ = ap.𝓒[1].z
    x₂ = ap.𝓒[2].x
    y₂ = ap.𝓒[2].y
    z₂ = ap.𝓒[2].z
    return ((x₂-x₁)^2+(y₂-y₁)^2+(z₂-z₁)^2)^0.5
end

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

function set𝝭!(ap::Element{:Seg2},x::Node)
    𝝭 = x[:𝝭]
    𝝭[1] = 0.5*(1.0-x.ξ)
    𝝭[2] = 0.5*(1.0+x.ξ)
end

function set∇𝝭!(ap::Element{:Seg2},x::Node)
    𝐿 = get𝐿(ap)
    ∂𝝭∂x = x[:∂𝝭∂x]
    ∂𝝭∂x[1] = -1.0/𝐿
    ∂𝝭∂x[2] = 1.0/𝐿
end