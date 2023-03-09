
function set𝝭!(ap::TRElement{:Tri3},x::Node)
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

function set∇𝝭!(ap::TRElement{:Tri3},x::Node)
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