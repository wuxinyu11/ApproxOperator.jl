
## Kirchhoff-Love plate
function (op::Operator{:∫κᵢⱼMᵢⱼdΩ})(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    D = op.D
    ν = op.ν
    for ξ in 𝓖
        _,_,_,B₁₁,B₁₂,B₂₂ = get∇²𝝭(ap,ξ)
        𝑤 = get𝑤(ap,ξ)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.id
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.id
                k[I,J] += D*(B₁₁[i]*B₁₁[j] + ν*(B₁₁[i]*B₂₂[j] + B₂₂[i]*B₁₁[j]) + B₂₂[i]*B₂₂[j] + 2*(1-ν)*B₁₂[i]*B₁₂[j])*𝑤
            end
        end
    end
end

function (op::Operator{:∫wqdΩ})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = get𝑤(ap,ξ)
        N = get𝝭(ap,ξ)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.id
            f[I] += N[i]*ξ.q*𝑤
        end
    end
end

function (op::Operator{:∫wVdΓ})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = get𝑤(ap,ξ)
        N = get𝝭(ap,ξ)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.id
            f[I] += N[i]*ξ.V*𝑤
        end
    end
end

function (op::Operator{:∫θₙMₙₙdΓ})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    n₁,n₂ = get𝒏(ap)
    for ξ in 𝓖
        𝑤 = get𝑤(ap,ξ)
        _,B₁,B₂ = get∇𝝭(ap,ξ)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.id
            f[I] -= (B₁[i]*n₁+B₂[i]*n₂)*ξ.M*𝑤
        end
    end
end

function (op::Operator{:∫∇wMdΓ})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    n₁,n₂ = get𝒏(ap)
    for ξ in 𝓖
        𝑤 = get𝑤(ap,ξ)
        _,B₁,B₂ = get∇𝝭(ap,ξ)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.id
            f[I] -= (B₁[i]*ξ.M₁+B₂[i]*ξ.M₂)*𝑤
        end
    end
end

function (op::Operator{:∫∇𝑛vθdΓ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    n₁,n₂ = get𝒏(ap)
    for ξ in 𝓖
        𝑤 = get𝑤(ap,ξ)
        _,B₁,B₂ = get∇𝝭(ap,ξ)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.id
            θᵢ = B₁[i]*n₁+B₂[i]*n₂
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.id
                θⱼ = B₁[j]*n₁+B₂[j]*n₂
                k[I,J] += op.α*θᵢ*θⱼ*𝑤
            end
            f[I] += op.α*θᵢ*ξ.θ*𝑤
        end
    end
end

function (op::Operator{:∫MₙₙθdΓ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    n₁,n₂,s₁,s₂ = get𝒏(ap)
    α = op.α
    D = op.D
    ν = op.ν
    D₁₁ = D*(n₁^2+ν*n₂^2)
    D₁₂ = 2*D*n₁*n₂*(1-ν)
    D₂₂ = D*(ν*n₁^2+n₂^2)
    for ξ in 𝓖
        𝑤 = get𝑤(ap,ξ)
        _,B₁,B₂,B₁₁,B₁₂,B₂₂ = get∇²𝝭(ap,ξ)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.id
            θᵢ = B₁[i]*n₁ + B₂[i]*n₂
            Mᵢ = D₁₁*B₁₁[i] + D₁₂*B₁₂[i] + D₂₂*B₂₂[i]
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.id
                θⱼ = B₁[j]*n₁ + B₂[j]*n₂
                Mⱼ = D₁₁*B₁₁[j] + D₁₂*B₁₂[j] + D₂₂*B₂₂[j]
                k[I,J] += (-Mᵢ*θⱼ-θᵢ*Mⱼ+α*θᵢ*θⱼ)*𝑤
            end
            f[I] += (-Mᵢ+α*θᵢ)*ξ.θ*𝑤
        end
    end
end

function (op::Operator{:∫VgdΓ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    n₁,n₂,s₁,s₂ = get𝒏(ap)
    α = op.α
    D = op.D
    ν = op.ν
    D₁₁₁ = -D*(n₁ + n₁*s₁*s₁ + ν*n₂*s₁*s₂)
    D₁₁₂ = -D*(n₂ + n₂*s₁*s₁ + 2*n₁*s₁*s₂ + (n₂*s₂*s₂ - n₂*s₁*s₁ - n₁*s₁*s₂)*ν)
    D₁₂₂ = -D*(n₁ + n₁*s₂*s₂ + 2*n₂*s₁*s₂ + (n₁*s₁*s₁ - n₁*s₂*s₂ - n₂*s₁*s₂)*ν)
    D₂₂₂ = -D*(n₂ + n₂*s₂*s₂ + ν*n₁*s₁*s₂)
    for ξ in 𝓖
        𝑤 = get𝑤(ap,ξ)
        N,_,_,_,_,_,B₁₁₁,B₁₁₂,B₁₂₂,B₂₂₂ = get∇³𝝭(ap,ξ)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.id
            Vᵢ = D₁₁₁*B₁₁₁[i] + D₁₁₂*B₁₁₂[i] + D₁₂₂*B₁₂₂[i] + D₂₂₂*B₂₂₂[i]
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.id
                Vⱼ = D₁₁₁*B₁₁₁[j] + D₁₁₂*B₁₁₂[j] + D₁₂₂*B₁₂₂[j] + D₂₂₂*B₂₂₂[j]
                k[I,J] += (-Vᵢ*N[j]-N[i]*Vⱼ+α*N[i]*N[j])*𝑤
            end
            f[I] += (-Vᵢ+α*N[i])*ξ.g*𝑤
        end
    end
end

function (op::Operator{:∫ṼgdΓ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    n₁,n₂,s₁,s₂ = get𝒏(ap)
    α = op.α
    D = op.D
    ν = op.ν
    D₁₁₁ = -D*(n₁ + n₁*s₁*s₁ + ν*n₂*s₁*s₂)
    D₁₁₂ = -D*(n₁*s₁*s₂ + (n₂*s₂*s₂ + n₂)*ν)
    D₁₂₁ = -D*(n₂ + n₂*s₁*s₁ + n₁*s₁*s₂ + (-n₂ - n₂*s₁*s₁ - n₁*s₁*s₂)*ν)
    D₁₂₂ = -D*(n₁ + n₁*s₂*s₂ + n₂*s₁*s₂ + (-n₁ - n₁*s₂*s₂ - n₂*s₁*s₂)*ν)
    D₂₂₁ = -D*(n₂*s₁*s₂ + (n₁*s₁*s₁ + n₁)*ν)
    D₂₂₂ = -D*(n₂ + n₂*s₂*s₂ + ν*n₁*s₁*s₂)
    for ξ in 𝓖
        𝑤 = get𝑤(ap,ξ)
        N,B₁₁₁,B₁₁₂,B₁₂₁,B₁₂₂,B₂₂₁,B₂₂₂ = get∇∇̃²𝝭(ap,ξ)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.id
            Vᵢ = D₁₁₁*B₁₁₁[i] + D₁₁₂*B₁₁₂[i] + D₁₂₁*B₁₂₁[i] + D₁₂₂*B₁₂₂[i] + D₂₂₁*B₂₂₁[i] + D₂₂₂*B₂₂₂[i]
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.id
                Vⱼ = D₁₁₁*B₁₁₁[j] + D₁₁₂*B₁₁₂[j] + D₁₂₁*B₁₂₁[j] + D₁₂₂*B₁₂₂[j] + D₂₂₁*B₂₂₁[j] + D₂₂₂*B₂₂₂[j]
                k[I,J] += (-Vᵢ*N[j]-N[i]*Vⱼ+α*N[i]*N[j])*𝑤
            end
            f[I] += (-Vᵢ+α*N[i])*ξ.g*𝑤
        end
    end
end

function (op::Operator{:wΔMₙₛ})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement{:Poi1}
    𝓒 = ap.𝓒; ξ = ap.𝓖[1]
    N = get𝝭(ap,ξ)
    for (i,xᵢ) in enumerate(𝓒)
        I = xᵢ.id
        f[I] -= N[i]*ξ.ΔM
    end
end

function (op::Operator{:ΔMₙₛg})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement{:Poi1}
    𝓒 = ap.𝓒; ξ = ap.𝓖[1]
    D = op.D
    ν = op.ν
    α = op.α
    Δn₁s₁ = ξ.Δn₁s₁
    Δn₁s₂n₂s₁ = ξ.Δn₁s₂n₂s₁
    Δn₂s₂ = ξ.Δn₂s₂
    D₁₁ = - D*(Δn₁s₁+Δn₂s₂*ν)
    D₁₂ = - D*(1-ν)*Δn₁s₂n₂s₁
    D₂₂ = - D*(Δn₁s₁*ν+Δn₂s₂)
    N,_,_,B₁₁,B₁₂,B₂₂ = get∇²𝝭(ap,ξ)
    for (i,xᵢ) in enumerate(𝓒)
        I = xᵢ.id
        ΔMₙₛᵢ = D₁₁*B₁₁[i] + D₁₂*B₁₂[i] + D₂₂*B₂₂[i]
        for (j,xⱼ) in enumerate(𝓒)
            J = xⱼ.id
            ΔMₙₛⱼ = D₁₁*B₁₁[j] + D₁₂*B₁₂[j] + D₂₂*B₂₂[j]
            k[I,J] += ΔMₙₛᵢ*N[j] + N[i]*ΔMₙₛⱼ + α*N[i]*N[j]
        end
        f[I] += (ΔMₙₛᵢ + α*N[i])*ξ.g
    end
end
