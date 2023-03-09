
@inline get𝒙(ap::T,ξ::Node) where T<:AbstractElement{:Tri3} = get𝒙(ap,ξ.ξ,ξ.η)

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

@inline get𝐽(ap::T,::Any) where T<:AbstractElement{:Tri3} = 2.0*get𝐴(ap)
@inline get𝑤(ap::T,ξ::Node) where T<:AbstractElement{:Tri3} = get𝐴(ap)*ξ.w

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

function set𝝭!(ap::Element{:Tri3},x::Node)
    𝝭 = x[:𝝭]
    𝝭[1] = x.ξ
    𝝭[2] = x.η
    𝝭[3] = 1.0-x.ξ-x.η
end
function set∇𝝭!(ap::Element{:Tri3},x::Node)
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