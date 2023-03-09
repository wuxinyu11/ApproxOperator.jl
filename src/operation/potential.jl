
function (op::Operator{:∫vudΩ})(ap::T,m::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ξ[:𝝭]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                m[I,J] += N[i]*N[j]*𝑤
            end
        end
    end
end

function (op::Operator{:∫∇v∇uvbdΩ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    kᶜ = op.k
    for ξ in 𝓖
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z]
        𝑤 = ξ.𝑤
        b = ξ.b
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += kᶜ*(B₁[i]*B₁[j] + B₂[i]*B₂[j] + B₃[i]*B₃[j])*𝑤
            end
            f[I] += N[i]*b*𝑤
        end
    end
end

function (op::Operator{:∫vₓuₓdx})(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    EA = op.EA
    for ξ in 𝓖
        B = ξ[:∂𝝭∂x]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += B[i]*EA*B[j]*𝑤
            end
        end
    end
end

function (op::Operator{:∫∇v∇udΩ})(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    kᶜ = op.k
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += kᶜ*(B₁[i]*B₁[j] + B₂[i]*B₂[j] + B₃[i]*B₃[j])*𝑤
            end
        end
    end
end

function (op::Operator{:∫∫∇v∇udxdy})(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    kᶜ = op.k
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += kᶜ*(B₁[i]*B₁[j] + B₂[i]*B₂[j])*𝑤
            end
        end
    end
end

function (op::Operator{:∫vbdΩ})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        b = ξ.b
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[I] += N[i]*b*𝑤
        end
    end
end

function (op::Operator{:∫vtdΓ})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        t = ξ.t
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            f[I] += N[i]*t*𝑤
        end
    end
end

function (op::Operator{:∫vgdΓ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    α = op.α
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += α*N[i]*N[j]*𝑤
            end
            f[I] += α*N[i]*g*𝑤
        end
    end
end

function (op::Operator{:∫λgdΓ})(ap1::T,ap2::S,g::AbstractMatrix{Float64},q::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    for j in 1:length(ap1.𝓖)
        ξ₁ = ap1.𝓖[j]
        ξ₂ = ap2.𝓖[j]
        𝑤 = ξ₁.𝑤
        N = ξ₁[:𝝭]
        N̄ = ξ₂[:𝝭]
        ḡ = ξ₁.g
        for (k,xₖ) in enumerate(ap2.𝓒)
            K = xₖ.𝐼
            for (i,xᵢ) in enumerate(ap1.𝓒)
                I = xᵢ.𝐼
                g[I,K] -= N[i]*N̄[k]*𝑤
            end
            q[K] -= N̄[k]*ḡ*𝑤
        end
    end
end

function (op::Operator{:∫λₙgdΓ})(ap1::T,ap2::S,g::AbstractMatrix{Float64},q::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    for j in 1:length(ap1.𝓖)
        ξ₁ = ap1.𝓖[j]
        ξ₂ = ap2.𝓖[j]
        𝑤 = ξ₁.𝑤
        N = ξ₁[:𝝭]
        N̄ = ξ₂[:𝝭]
        ḡ = ξ₁.g
        sn = sign(ξ₁.n₁ + ξ₂.n₂)
        for (k,xₖ) in enumerate(ap2.𝓒)
            K = xₖ.𝐼
            for (i,xᵢ) in enumerate(ap1.𝓒)
                I = xᵢ.𝐼
                g[I,K] -= sn*N[i]*N̄[k]*𝑤
            end
            q[K] -= sn*N̄[k]*ḡ*𝑤
        end
    end
end

function (op::Operator{:∫∇𝑛vgdΓ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    kᶜ = op.k
    α = op.α
    for ξ in 𝓖
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        B₃ = ξ[:∂𝝭∂z]
        𝑤 = ξ.𝑤
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        n₃ = ξ.n₃
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (-kᶜ*((B₁[i]*n₁+B₂[i]*n₂+B₃[i]*n₃)*N[j]+N[i]*(B₁[j]*n₁+B₂[j]*n₂+B₃[j]*n₃)) + α*N[i]*N[j])*𝑤
            end
            f[I] += (-kᶜ*(B₁[i]*n₁+B₂[i]*n₂+B₃[i]*n₃) + α*N[i])*g*𝑤
        end
    end
end

function (op::Operator{:∫∇𝑛vgds})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒;𝓖 = ap.𝓖
    kᶜ = op.k
    α = op.α
    for ξ in 𝓖
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x]
        B₂ = ξ[:∂𝝭∂y]
        𝑤 = ξ.𝑤
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (-kᶜ*((B₁[i]*n₁+B₂[i]*n₂)*N[j]+N[i]*(B₁[j]*n₁+B₂[j]*n₂)) + α*N[i]*N[j])*𝑤
            end
            f[I] += (-kᶜ*(B₁[i]*n₁+B₂[i]*n₂) + α*N[i])*g*𝑤
        end
    end
end

function (op::Operator{:∫∇̄𝑛vgdΓ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    kᶜ = op.k
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        B₁ = ξ[:∂𝝭∂x_]
        B₂ = ξ[:∂𝝭∂y_]
        B₃ = ξ[:∂𝝭∂z_]
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        n₃ = ξ.n₃
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += kᶜ*(B₁[i]*n₁+B₂[i]*n₂+B₃[i]*n₃)*N[j]*𝑤
            end
            f[I] += kᶜ*(B₁[i]*n₁+B₂[i]*n₂+B₃[i]*n₃)*g*𝑤
        end
    end
end

function (op::Operator{:g})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64};dof::Symbol=:d) where T<:AbstractElement{:Poi1}
    x = ap.𝓒[1]
    j = x.𝐼
    g = getproperty(x,dof)
    for i in 1:length(f)
        f[i] -= k[i,j]*g
    end
    k[j,:] .= 0.
    k[:,j] .= 0.
    k[j,j] = 1.
    f[j] = g
end
