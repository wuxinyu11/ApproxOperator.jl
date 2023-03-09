
function (op::Operator{:∫κᵢⱼMᵢⱼdΩ})(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    D = op.D
    ν = op.ν
    for ξ in 𝓖
        B₁₁ = ξ[:∂²𝝭∂x²]
        B₁₂ = ξ[:∂²𝝭∂x∂y]
        B₂₂ = ξ[:∂²𝝭∂y²]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += D*(B₁₁[i]*B₁₁[j] + ν*(B₁₁[i]*B₂₂[j] + B₂₂[i]*B₁₁[j]) + B₂₂[i]*B₂₂[j] + 2*(1-ν)*B₁₂[i]*B₁₂[j])*𝑤
            end
        end
    end
end

function (op::Operator{:∫κ̃ᵢⱼM̃ᵢⱼdΩ})(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    D = op.D
    ν = op.ν
    for ξ in 𝓖
        B₁₁ = ξ[:∂²𝝭∂x²]
        B₁₂ = ξ[:∂²𝝭∂x∂y]
        B₂₁ = ξ[:∂²𝝭∂y∂x]
        B₂₂ = ξ[:∂²𝝭∂y²]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += D*(B₁₁[i]*B₁₁[j] + ν*(B₁₁[i]*B₂₂[j] + B₂₂[i]*B₁₁[j]) + B₂₂[i]*B₂₂[j] + (1-ν)*(B₁₂[i]*B₁₂[j]+B₂₁[i]*B₂₁[j]))*𝑤
            end
        end
    end
end

function (op::Operator{:∫wqdΩ})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        q = ξ.q
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[I] += N[i]*q*𝑤
        end
    end
end

function (op::Operator{:∫wVdΓ})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        V = ξ.V
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[I] += N[i]*V*𝑤
        end
    end
end

function (op::Operator{:∫θₙMₙₙdΓ})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        M = ξ.M
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[I] -= (B₁[i]*n₁+B₂[i]*n₂)*M*𝑤
        end
    end
end

function (op::Operator{:∫∇wMdΓ})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    n₁,n₂ = get𝒏(ap)
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        M₁ = ξ.M₁
        M₂ = ξ.M₂
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[I] -= (B₁[i]*M₁+B₂[i]*M₂)*𝑤
        end
    end
end

function (op::Operator{:∫∇𝑛vθdΓ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    α = op.α
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        θ = ξ.θ
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            θᵢ = B₁[i]*n₁+B₂[i]*n₂
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                θⱼ = B₁[j]*n₁+B₂[j]*n₂
                k[I,J] += α*θᵢ*θⱼ*𝑤
            end
            f[I] += α*θᵢ*θ*𝑤
        end
    end
end

function (op::Operator{:∫MₙₙθdΓ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    α = op.α
    D = op.D
    ν = op.ν
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        D₁₁ = -D*(n₁^2+ν*n₂^2)
        D₁₂ = -2*D*n₁*n₂*(1-ν)
        D₂₂ = -D*(ν*n₁^2+n₂^2)
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₁₁ = ξ[:∂²𝝭∂x²]
        B₁₂ = ξ[:∂²𝝭∂x∂y]
        B₂₂ = ξ[:∂²𝝭∂y²]
        θ = ξ.θ
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            θᵢ = B₁[i]*n₁ + B₂[i]*n₂
            Mᵢ = D₁₁*B₁₁[i] + D₁₂*B₁₂[i] + D₂₂*B₂₂[i]
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                θⱼ = B₁[j]*n₁ + B₂[j]*n₂
                Mⱼ = D₁₁*B₁₁[j] + D₁₂*B₁₂[j] + D₂₂*B₂₂[j]
                k[I,J] += (Mᵢ*θⱼ+θᵢ*Mⱼ+α*θᵢ*θⱼ)*𝑤
            end
            f[I] += (Mᵢ+α*θᵢ)*θ*𝑤
        end
    end
end

function (op::Operator{:∫VgdΓ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    α = op.α
    D = op.D
    ν = op.ν
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        s₁ = -ξ.n₂
        s₂ = ξ.n₁
        D₁₁₁ = -D*(n₁ + n₁*s₁*s₁ + ν*n₂*s₁*s₂)
        D₁₁₂ = -D*(n₂ + n₂*s₁*s₁ + 2*n₁*s₁*s₂ + (n₂*s₂*s₂ - n₂*s₁*s₁ - n₁*s₁*s₂)*ν)
        D₁₂₂ = -D*(n₁ + n₁*s₂*s₂ + 2*n₂*s₁*s₂ + (n₁*s₁*s₁ - n₁*s₂*s₂ - n₂*s₁*s₂)*ν)
        D₂₂₂ = -D*(n₂ + n₂*s₂*s₂ + ν*n₁*s₁*s₂)
        N = ξ[:𝝭]
        B₁₁₁ = ξ[:∂³𝝭∂x³]
        B₁₁₂ = ξ[:∂³𝝭∂x²∂y]
        B₁₂₂ = ξ[:∂³𝝭∂x∂y²]
        B₂₂₂ = ξ[:∂³𝝭∂y³]
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            Vᵢ = D₁₁₁*B₁₁₁[i] + D₁₁₂*B₁₁₂[i] + D₁₂₂*B₁₂₂[i] + D₂₂₂*B₂₂₂[i]
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                Vⱼ = D₁₁₁*B₁₁₁[j] + D₁₁₂*B₁₁₂[j] + D₁₂₂*B₁₂₂[j] + D₂₂₂*B₂₂₂[j]
                k[I,J] += (-Vᵢ*N[j]-N[i]*Vⱼ+α*N[i]*N[j])*𝑤
            end
            f[I] += (-Vᵢ+α*N[i])*g*𝑤
        end
    end
end

function (op::Operator{:∫M̃ₙₙθdΓ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    n₁ = 𝓖[1].n₁
    n₂ = 𝓖[1].n₂
    s₁ = 𝓖[1].s₁
    s₂ = 𝓖[1].s₂
    D = op.D
    ν = op.ν
    D₁₁ = -D*(n₁^2+ν*n₂^2)
    D₁₂ = -2*D*n₁*n₂*(1-ν)
    D₂₂ = -D*(ν*n₁^2+n₂^2)
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₁₁ = ξ[:∂²𝝭∂x²]
        B₁₂ = ξ[:∂²𝝭∂x∂y]
        B₂₂ = ξ[:∂²𝝭∂y²]
        B̄₁₁ = ξ[:∂²𝝭∂x²_]
        B̄₁₂ = ξ[:∂²𝝭∂x∂y_]
        B̄₂₂ = ξ[:∂²𝝭∂y²_]
        θ = ξ.θ
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            θᵢ = B₁[i]*n₁ + B₂[i]*n₂
            Mᵢ = D₁₁*B₁₁[i] + D₁₂*B₁₂[i] + D₂₂*B₂₂[i]
            M̄ᵢ = D₁₁*B̄₁₁[i] + D₁₂*B̄₁₂[i] + D₂₂*B̄₂₂[i]
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                θⱼ = B₁[j]*n₁ + B₂[j]*n₂
                Mⱼ = D₁₁*B₁₁[j] + D₁₂*B₁₂[j] + D₂₂*B₂₂[j]
                k[I,J] += (Mᵢ*θⱼ+θᵢ*Mⱼ-M̄ᵢ*θⱼ)*𝑤
            end
            f[I] += (Mᵢ-M̄ᵢ)*θ*𝑤
        end
    end
end

function (op::Operator{:∫ṼgdΓ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    n₁ = 𝓖[1].n₁
    n₂ = 𝓖[1].n₂
    s₁ = 𝓖[1].s₁
    s₂ = 𝓖[1].s₂
    D = op.D
    ν = op.ν
    D₁₁₁ = -D*(n₁ + n₁*s₁*s₁ + ν*n₂*s₁*s₂)
    D₁₁₂ = -D*(n₁*s₁*s₂ + (n₂*s₂*s₂ + n₂)*ν)
    D₁₂₁ = -D*(n₂ + n₂*s₁*s₁ + n₁*s₁*s₂ + (-n₂ - n₂*s₁*s₁ - n₁*s₁*s₂)*ν)
    D₁₂₂ = -D*(n₁ + n₁*s₂*s₂ + n₂*s₁*s₂ + (-n₁ - n₁*s₂*s₂ - n₂*s₁*s₂)*ν)
    D₂₂₁ = -D*(n₂*s₁*s₂ + (n₁*s₁*s₁ + n₁)*ν)
    D₂₂₂ = -D*(n₂ + n₂*s₂*s₂ + ν*n₁*s₁*s₂)
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁₁₁ = ξ[:∂∂²𝝭∂x²∂x]
        B₁₁₂ = ξ[:∂∂²𝝭∂x²∂y]
        B₁₂₁ = ξ[:∂∂²𝝭∂x∂y∂x]
        B₁₂₂ = ξ[:∂∂²𝝭∂x∂y∂y]
        B₂₂₁ = ξ[:∂∂²𝝭∂y²∂x]
        B₂₂₂ = ξ[:∂∂²𝝭∂y²∂y]
        B̄₁₁₁ = ξ[:∂∂²𝝭∂x²∂x_]
        B̄₁₁₂ = ξ[:∂∂²𝝭∂x²∂y_]
        B̄₁₂₁ = ξ[:∂∂²𝝭∂x∂y∂x_]
        B̄₁₂₂ = ξ[:∂∂²𝝭∂x∂y∂y_]
        B̄₂₂₁ = ξ[:∂∂²𝝭∂y²∂x_]
        B̄₂₂₂ = ξ[:∂∂²𝝭∂y²∂y_]
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            Vᵢ = D₁₁₁*B₁₁₁[i] + D₁₁₂*B₁₁₂[i] + D₁₂₁*B₁₂₁[i] + D₁₂₂*B₁₂₂[i] + D₂₂₁*B₂₂₁[i] + D₂₂₂*B₂₂₂[i]
            V̄ᵢ = D₁₁₁*B̄₁₁₁[i] + D₁₁₂*B̄₁₁₂[i] + D₁₂₁*B̄₁₂₁[i] + D₁₂₂*B̄₁₂₂[i] + D₂₂₁*B̄₂₂₁[i] + D₂₂₂*B̄₂₂₂[i]
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                Vⱼ = D₁₁₁*B₁₁₁[j] + D₁₁₂*B₁₁₂[j] + D₁₂₁*B₁₂₁[j] + D₁₂₂*B₁₂₂[j] + D₂₂₁*B₂₂₁[j] + D₂₂₂*B₂₂₂[j]
                k[I,J] -= (Vᵢ*N[j]+N[i]*Vⱼ-V̄ᵢ*N[j])*𝑤
            end
            f[I] -= (Vᵢ-V̄ᵢ)*g*𝑤
        end
    end
end

function (op::Operator{:wΔMₙₛ})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        ΔM = ξ.ΔM
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[I] -= N[i]*ΔM
        end
    end
end

function (op::Operator{:ΔMₙₛg})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    D = op.D
    ν = op.ν
    α = op.α
    for ξ in 𝓖
        Δn₁s₁ = ξ.Δn₁s₁
        Δn₁s₂n₂s₁ = ξ.Δn₁s₂n₂s₁
        Δn₂s₂ = ξ.Δn₂s₂
        D₁₁ = - D*(Δn₁s₁+Δn₂s₂*ν)
        D₁₂ = - D*(1-ν)*Δn₁s₂n₂s₁
        D₂₂ = - D*(Δn₁s₁*ν+Δn₂s₂)
        N = ξ[:𝝭]
        B₁₁ = ξ[:∂²𝝭∂x²]
        B₁₂ = ξ[:∂²𝝭∂x∂y]
        B₂₂ = ξ[:∂²𝝭∂y²]
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            ΔMₙₛᵢ = D₁₁*B₁₁[i] + D₁₂*B₁₂[i] + D₂₂*B₂₂[i]
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                ΔMₙₛⱼ = D₁₁*B₁₁[j] + D₁₂*B₁₂[j] + D₂₂*B₂₂[j]
                k[I,J] += ΔMₙₛᵢ*N[j] + N[i]*ΔMₙₛⱼ + α*N[i]*N[j]
            end
            f[I] += (ΔMₙₛᵢ + α*N[i])*g
        end
    end
end

function (op::Operator{:ΔM̃ₙₛg})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    D = op.D
    ν = op.ν
    for ξ in 𝓖
        Δn₁s₁ = ξ.Δn₁s₁
        Δn₁s₂n₂s₁ = ξ.Δn₁s₂n₂s₁
        Δn₂s₂ = ξ.Δn₂s₂
        D₁₁ = - D*(Δn₁s₁+Δn₂s₂*ν)
        D₁₂ = - D*(1-ν)*Δn₁s₂n₂s₁
        D₂₂ = - D*(Δn₁s₁*ν+Δn₂s₂)
        N = ξ[:𝝭]
        B₁₁ = ξ[:∂²𝝭∂x²]
        B₁₂ = ξ[:∂²𝝭∂x∂y]
        B₂₂ = ξ[:∂²𝝭∂y²]
        B̄₁₁ = ξ[:∂²𝝭∂x²_]
        B̄₁₂ = ξ[:∂²𝝭∂x∂y_]
        B̄₂₂ = ξ[:∂²𝝭∂y²_]
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            ΔMₙₛᵢ = D₁₁*B₁₁[i] + D₁₂*B₁₂[i] + D₂₂*B₂₂[i]
            ΔM̄ₙₛᵢ = D₁₁*B̄₁₁[i] + D₁₂*B̄₁₂[i] + D₂₂*B̄₂₂[i]
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                ΔMₙₛⱼ = D₁₁*B₁₁[j] + D₁₂*B₁₂[j] + D₂₂*B₂₂[j]
                k[I,J] += ΔMₙₛᵢ*N[j] + N[i]*ΔMₙₛⱼ - ΔM̄ₙₛᵢ*N[j]
            end
            f[I] += (ΔMₙₛᵢ - ΔM̄ₙₛᵢ)*g
        end
    end
end
