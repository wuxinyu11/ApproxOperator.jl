
function checkIC(a::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Tri3},b::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Tri3}) where {𝒑,𝑠,𝜙}
    x₁ = a.𝓒[1].x;y₁ = a.𝓒[1].y
    x₂ = a.𝓒[2].x;y₂ = a.𝓒[2].y
    x₃ = a.𝓒[3].x;y₃ = a.𝓒[3].y
    n₁₁ = y₃-y₂;n₂₁ = y₁-y₃;n₃₁ = y₂-y₁
    n₁₂ = x₂-x₃;n₂₂ = x₃-x₁;n₃₂ = x₁-x₂
    nᵈ = length(a.𝓒)
    nᵖ = get𝑛𝒒(a)
    # nᵖ = get𝑛𝒑(a)
    f₁ = zeros(nᵈ,nᵖ)
    f₂ = zeros(nᵈ,nᵖ)
    for ξ in a.𝓖
        N,B₁,B₂ = get∇𝝭(a,ξ)
        𝒙 = get𝒙(a,ξ)
        𝒒 = get𝒒(a,ξ)
        p = get𝒑(a,𝒙)
        𝑤 = get𝑤(a,ξ)
        for i in 1:nᵈ
            for j in 1:nᵖ
                f₁[i,j] += B₁[i]*p[j]*𝑤
                f₂[i,j] += B₂[i]*p[j]*𝑤
            end
        end
    end
    for ξ in b.𝓖
        N = get𝝭(b,ξ)
        𝒒,∂𝒒∂ξ,∂𝒒∂η = get∇𝒒(b,ξ)
        𝒙 = get𝒙(b,ξ)
        p,∂p∂x,∂p∂y = get∇𝒑(b,𝒙)
        wᵇ = ξ.wᵇ
        w = ξ.w
        𝑤 = get𝑤(b,ξ)
        nᵇ₁ = 0.0;nᵇ₂ = 0.0
        ξ.ξ == 0.0 ? (nᵇ₁ += n₁₁;nᵇ₂ += n₁₂) : nothing
        ξ.η == 0.0 ? (nᵇ₁ += n₂₁;nᵇ₂ += n₂₂) : nothing
        ξ.ξ+ξ.η ≈ 1.0 ? (nᵇ₁ += n₃₁;nᵇ₂ += n₃₂) : nothing
        for j in 1:nᵖ
            W₁ = p[j]*nᵇ₁*wᵇ - ∂p∂x[j]*𝑤
            W₂ = p[j]*nᵇ₂*wᵇ - ∂p∂y[j]*𝑤
            for i in 1:nᵈ
                f₁[i,j] -= N[i]*W₁
                f₂[i,j] -= N[i]*W₂
            end
        end
    end
    return f₁,f₂
end

function checkCC(a::ReproducingKernel)
    nᵖ = get𝑛𝒑(a)
    nⁱ = length(a.𝓖)
    f = zeros(nⁱ,nᵖ)
    f₁ = zeros(nⁱ,nᵖ)
    f₂ = zeros(nⁱ,nᵖ)
    for i in 1:length(a.𝓖)
        ξ = a.𝓖[i]
        𝒙 = get𝒙(a,ξ)
        p,∂p∂x,∂p∂y = get∇𝒑(a,𝒙)
        N,B₁,B₂ = get∇𝝭(a,ξ)
        for j in 1:nᵖ
            for k in 1:length(a.𝓒)
                𝒙̄ = (a.𝓒[k].x,a.𝓒[k].y,a.𝓒[k].z)
                p̄ = get𝒑(a,𝒙̄)
                f[i,j]  += N[k]*p̄[j]
                f₁[i,j] += B₁[k]*p̄[j]
                f₂[i,j] += B₂[k]*p̄[j]
            end
            f[i,j]  -= p[j]
            f₁[i,j] -= ∂p∂x[j]
            f₂[i,j] -= ∂p∂y[j]
        end
    end
    return f,f₁,f₂
end

function checkConsistency(a::ReproducingKernel,get𝝭::Function,get𝒑::Function)
    nᵖ = get𝑛𝒑(a)
    nⁱ = length(a.𝓖)
    f = zeros(nᵖ)
    𝒑ʰ = zeros(nᵖ)
    for ξ in (a.𝓖)
        𝒙 = get𝒙(a,ξ)
        𝝭 = get𝝭(a,ξ)
        𝑤 = get𝑤(a,ξ)
        𝒑 = get𝒑(a,𝒙)
        fill!(𝒑,0.0)
        for (j,node) in enumerate(a.𝓒)
            𝒙ᵢ = (node.x,node.y,node.z)
            𝒑ᵢ = get𝒑(a,𝒙ᵢ)
            for i = 1:nᵖ
                𝒑ʰ[i] += 𝝭[i][j]*𝒑ᵢ[j]
            end
        end
        f .+= (𝒑 .- 𝒑ʰ).^2 .* 𝑤
    end
    return f
end
