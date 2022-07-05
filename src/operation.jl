
## Counstruction
struct Operator{T,D}
    type::Val{T}
    data::Dict{Symbol,D}
end
Operator(t::Symbol) = Operator(Val(t),Dict{Symbol,Any}())
function Operator(t::Symbol,d::Pair{Symbol,D}...) where D<:Any
    return Operator(Val(t),Dict(d))
end


## General Functions
push!(op::Operator,d::Pair{Symbol,D}...) where D<:Any = push!(op.data,d...)
@inline getproperty(op::Operator,f::Symbol) = hasfield(Operator,f) ? getfield(op,f) : getfield(op,:data)[f]
@inline function setproperty!(op::Operator,f::Symbol,x)
    getfield(op,:data)[f] = x
end

@inline function (op::Operator)(aps::Vector{T},gps::Vector{S},k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    for i in 1:length(aps)
        @inbounds op(aps[i],gps[i],k,f)
    end
end

@inline function (op::Operator)(aps::Vector{T},k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    for ap in aps
        op(ap,k,f)
    end
end
@inline function (op::Operator)(aps::Vector{T},k::AbstractMatrix{Float64}) where T<:AbstractElement
    for ap in aps
        op(ap,k)
    end
end
@inline function (op::Operator)(aps::Vector{T},f::AbstractVector{Float64}) where T<:AbstractElement
    for ap in aps
        op(ap,f)
    end
end

@inline function (op::Operator)(aps::Vector{T},s::Symbol) where T<:AbstractElement
    for ap in aps
        op(ap,s)
    end
end
@inline function (op::Operator)(aps::Vector{T}) where T<:AbstractElement
    for ap in aps
        op(ap)
    end
end

function prescribe!(ap::T,s::Symbol,f::Function) where T<:AbstractElement
    𝓖 = ap.𝓖
    for ξ in 𝓖
        𝒙 = get𝒙(ap,ξ)
        if applicable(f,𝒙...)
            v = f(𝒙...)
        else
            𝒏 = get𝒏(ap)
            v = f(𝒙...,𝒏...)
        end
        setproperty!(ξ,s,v)
    end
end

function prescribe!(aps::Vector{T},s::Symbol,f::Function) where T<:AbstractElement
    for ap in aps
        prescribe!(ap,s,f)
    end
end

##
function (op::Operator{:𝑓𝑣})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ.𝝭
        u = ξ.u
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.I
            f[I] += N[i]*u*𝑤
        end
    end
end

## Potential Problem
function (op::Operator{:∫∇v∇uvbdΩ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N,B₁,B₂,B₃ = get∇𝝭(ap,ξ)
        𝑤 = get𝑤(ap,ξ)
        for i in 1:length(𝓒)
            I = 𝓒[i].id
            for j in 1:length(𝓒)
                J = 𝓒[j].id
                k[I,J] += op.k*(B₁[i]*B₁[j] + B₂[i]*B₂[j] + B₃[i]*B₃[j])*𝑤
            end
            f[I] += N[i]*ξ.b*𝑤
        end
    end
end

function (op::Operator{:∫∇v∇udΩ})(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        _,B₁,B₂,B₃ = get∇𝝭(ap,ξ)
        𝑤 = get𝑤(ap,ξ)
        for i in 1:length(𝓒)
            I = 𝓒[i].id
            for j in 1:length(𝓒)
                J = 𝓒[j].id
                k[I,J] += op.k*(B₁[i]*B₁[j] + B₂[i]*B₂[j] + B₃[i]*B₃[j])*𝑤
            end
        end
    end
end

function (op::Operator{:∫vbdΩ})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = get𝑤(ap,ξ)
        N = get𝝭(ap,ξ)
        for i in 1:length(𝓒)
            I = 𝓒[i].id
            f[I] += N[i]*ξ.b*𝑤
        end
    end
end

function (op::Operator{:∫vtdΓ})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = get𝑤(ap,ξ)
        N = get𝝭(ap,ξ)
        for i in 1:length(𝓒)
            I = 𝓒[i].id
            f[I] += N[i]*ξ.t*𝑤
        end
    end
end

function (op::Operator{:∫vgdΓ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = get𝑤(ap,ξ)
        N = get𝝭(ap,ξ)
        for i in 1:length(𝓒)
            I = 𝓒[i].id
            for j in 1:length(𝓒)
                J = 𝓒[j].id
                k[I,J] += op.α*N[i]*N[j]*𝑤
            end
            f[I] += op.α*N[i]*ξ.g*𝑤
        end
    end
end

function (op::Operator{:∫λgdΓ})(ap1::T,ap2::S,g::AbstractMatrix{Float64},q::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    for ξ in ap1.𝓖
        𝑤 = get𝑤(ap1,ξ)
        N = get𝝭(ap1,ξ)
        N̄ = get𝝭(ap2,ξ)
        for k in 1:length(ap2.𝓒)
            K = ap2.𝓒[k].id
            for i in 1:length(ap1.𝓒)
                I = ap1.𝓒[i].id
                g[I,K] -= N[i]*N̄[k]*𝑤
            end
            q[K] -= N̄[k]*ξ.g*𝑤
        end
    end
end

function (op::Operator{:∫∇𝑛vgdΓ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    for ξ in ap.𝓖
        N,B₁,B₂,B₃ = get∇𝝭(ap,ξ)
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        n₃ = ξ.n₃
        for i in 1:length(ap.𝓒)
            I = ap.𝓒[i].id
            for j in 1:length(ap.𝓒)
                J = ap.𝓒[j].id
                k[I,J] -= op.k*((B₁[i]*n₁+B₂[i]*n₂+B₃[i]*n₃)*N[j] + N[i]*(B₁[j]*n₁+B₂[j]*n₂+B₃[j]*n₃))*ξ.w
            end
            f[I] -= op.k*(B₁[i]*n₁+B₂[i]*n₂+B₃[i]*n₃)*ξ.g*ξ.w
        end
    end
end

function (op::Operator{:∫∇̄𝑛vgdΓ})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = get𝑤(ap,ξ)
        N,B₁,B₂,B₃ = get∇̄𝝭(ap,ξ)
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        n₃ = ξ.n₃
        g = ξ.g
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.id
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.id
                k[I,J] += (B₁[i]*n₁+B₂[i]*n₂+B₃[i]*n₃)*N[j]*𝑤
            end
            f[I] += (B₁[i]*n₁+B₂[i]*n₂+B₃[i]*n₃)*ξ.g*𝑤
        end
    end
end

function (op::Operator{:g})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64};dof::Symbol=:d) where T<:AbstractElement{:Poi1}
    x = ap.𝓒[1]
    j = x.id
    g = getproperty(x,dof)
    for i in 1:length(f)
        f[i] -= k[i,j]*g
    end
    k[j,:] .= 0.
    k[:,j] .= 0.
    k[j,j] = 1.
    f[j] = g
end

## Plane Strain
function (op::Operator{:∫∫εᵢⱼσᵢⱼvᵢbᵢdxdy})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N,B₁,B₂ = get∇𝝭(ap,ξ)
        𝑤 = get𝑤(ap,ξ)
        E = op.E
        ν = op.ν
        Cᵢᵢᵢᵢ = E/(1-ν^2)
        Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
        for i in 1:length(𝓒)
            I = 𝓒[i].id
            for j in 1:length(𝓒)
                J = 𝓒[j].id
                k[2*I-1,2*J-1] += (Cᵢᵢᵢᵢ*B₁[i]*B₁[j] + Cᵢⱼᵢⱼ*B₂[i]*B₂[j])*𝑤
                k[2*I-1,2*J]   += (Cᵢᵢⱼⱼ*B₁[i]*B₂[j] + Cᵢⱼᵢⱼ*B₂[i]*B₁[j])*𝑤
                k[2*I,2*J-1]   += (Cᵢᵢⱼⱼ*B₂[i]*B₁[j] + Cᵢⱼᵢⱼ*B₁[i]*B₂[j])*𝑤
                k[2*I,2*J]     += (Cᵢᵢᵢᵢ*B₂[i]*B₂[j] + Cᵢⱼᵢⱼ*B₁[i]*B₁[j])*𝑤
            end
            f[2*I-1] += N[i]*ξ.b₁*𝑤
            f[2*I]   += N[i]*ξ.b₂*𝑤
        end
    end
end

function (op::Operator{:∫∫εᵢⱼσᵢⱼdxdy})(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        _,B₁,B₂ = get∇𝝭(ap,ξ)
        𝑤 = get𝑤(ap,ξ)
        E = op.E
        ν = op.ν
        Cᵢᵢᵢᵢ = E/(1-ν^2)
        Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
        for i in 1:length(𝓒)
            I = 𝓒[i].id
            for j in 1:length(𝓒)
                J = 𝓒[j].id
                k[2*I-1,2*J-1] += (Cᵢᵢᵢᵢ*B₁[i]*B₁[j] + Cᵢⱼᵢⱼ*B₂[i]*B₂[j])*𝑤
                k[2*I-1,2*J]   += (Cᵢᵢⱼⱼ*B₁[i]*B₂[j] + Cᵢⱼᵢⱼ*B₂[i]*B₁[j])*𝑤
                k[2*I,2*J-1]   += (Cᵢᵢⱼⱼ*B₂[i]*B₁[j] + Cᵢⱼᵢⱼ*B₁[i]*B₂[j])*𝑤
                k[2*I,2*J]     += (Cᵢᵢᵢᵢ*B₂[i]*B₂[j] + Cᵢⱼᵢⱼ*B₁[i]*B₁[j])*𝑤
            end
        end
    end
end

function (op::Operator{:∫∫εᵛᵢⱼσᵛᵢⱼdxdy})(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        _,B₁,B₂ = get∇𝝭(ap,ξ)
        𝑤 = get𝑤(ap,ξ)
        E = op.E
        ν = op.ν
        Cᵛ = E/(1-2*ν)
        for i in 1:length(𝓒)
            I = 𝓒[i].id
            for j in 1:length(𝓒)
                J = 𝓒[j].id
                k[2*I-1,2*J-1] += Cᵛ/3*B₁[i]*B₁[j]*𝑤
                k[2*I-1,2*J]   += Cᵛ/3*B₁[i]*B₂[j]*𝑤
                k[2*I,2*J-1]   += Cᵛ/3*B₂[i]*B₁[j]*𝑤
                k[2*I,2*J]     += Cᵛ/3*B₂[i]*B₂[j]*𝑤
            end
        end
    end
end

function (op::Operator{:∫∫εᵈᵢⱼσᵈᵢⱼdxdy})(ap::T,k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        _,B₁,B₂ = get∇𝝭(ap,ξ)
        𝑤 = get𝑤(ap,ξ)
        E = op.E
        ν = op.ν
        Cᵈ = E/(1+ν)
        for i in 1:length(𝓒)
            I = 𝓒[i].id
            for j in 1:length(𝓒)
                J = 𝓒[j].id
                k[2*I-1,2*J-1] += Cᵈ*( 2/3*B₁[i]*B₁[j]+1/2*B₂[i]*B₂[j])*𝑤
                k[2*I-1,2*J]   += Cᵈ*(-1/3*B₁[i]*B₂[j]+1/2*B₂[i]*B₁[j])*𝑤
                k[2*I,2*J-1]   += Cᵈ*(-1/3*B₂[i]*B₁[j]+1/2*B₁[i]*B₂[j])*𝑤
                k[2*I,2*J]     += Cᵈ*( 2/3*B₂[i]*B₂[j]+1/2*B₁[i]*B₁[j])*𝑤
            end
        end
    end
end

function (op::Operator{:∫∫vᵢbᵢdxdy})(ap::T,f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = get𝝭(ap,ξ)
        𝑤 = get𝑤(ap,ξ)
        for i in 1:length(𝓒)
            I = 𝓒[i].id
            f[2*I-1] += N[i]*ξ.b₁*𝑤
            f[2*I]   += N[i]*ξ.b₂*𝑤
        end
    end
end

function (op::Operator{:∫vᵢtᵢds})(ap::T,f::Vector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = get𝝭(ap,ξ)
        𝑤 = get𝑤(ap,ξ)
        for i in 1:length(𝓒)
            I = 𝓒[i].id
            f[2*I-1] += N[i]*ξ.t₁*𝑤
            f[2*I]   += N[i]*ξ.t₂*𝑤
        end
    end
end

function (op::Operator{:∫λᵢgᵢds})(ap1::T,ap2::S,g::AbstractMatrix{Float64},q::AbstractVector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    for ξ in ap1.𝓖
        𝑤 = get𝑤(ap1,ξ)
        N = get𝝭(ap1,ξ)
        N̄ = get𝝭(ap2,ξ)
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        n₁₁ = ξ.n₁₁
        n₂₂ = ξ.n₂₂
        n₁₂ = ξ.n₁₂
        for k in 1:length(ap2.𝓒)
            K = ap2.𝓒[k].id
            for i in 1:length(ap1.𝓒)
                I = ap1.𝓒[i].id
                N̄ₖNᵢ = N̄[k]*N[i]
                g[2*K-1,2*I-1] -= n₁₁*N̄ₖNᵢ*𝑤
                g[2*K-1,2*I]   -= n₁₂*N̄ₖNᵢ*𝑤
                g[2*K,2*I-1]   -= n₁₂*N̄ₖNᵢ*𝑤
                g[2*K,2*I]     -= n₂₂*N̄ₖNᵢ*𝑤
            end
            q[2*K-1] -= N̄[k]*(n₁₁*g₁+n₁₂*g₂)*𝑤
            q[2*K]   -= N̄[k]*(n₁₂*g₁+n₂₂*g₂)*𝑤
        end
    end
end

function (op::Operator{:∫σᵢⱼnⱼgᵢds})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N,B₁,B₂ = get∇𝝭(ap,ξ)
        E = op.E
        ν = op.ν
        Cᵢᵢᵢᵢ = E/(1-ν^2)
        Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
        n₁₁ = ξ.n₁₁
        n₁₂ = ξ.n₁₂
        n₂₂ = ξ.n₂₂
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        C₁₁₁ = Cᵢᵢᵢᵢ*n₁*n₁₁+Cᵢᵢⱼⱼ*n₂*n₁₂
        C₁₁₂ = Cᵢᵢᵢᵢ*n₁*n₁₂+Cᵢᵢⱼⱼ*n₂*n₂₂
        C₂₂₁ = Cᵢᵢⱼⱼ*n₁*n₁₁+Cᵢᵢᵢᵢ*n₂*n₁₂
        C₂₂₂ = Cᵢᵢⱼⱼ*n₁*n₁₂+Cᵢᵢᵢᵢ*n₂*n₂₂
        C₁₂₁ = Cᵢⱼᵢⱼ*(n₁*n₁₂+n₂*n₁₁)
        C₁₂₂ = Cᵢⱼᵢⱼ*(n₂*n₁₂+n₁*n₂₂)
        for i in 1:length(𝓒)
            I = 𝓒[i].id
            for j in 1:length(𝓒)
                J = 𝓒[j].id
                k[2*I-1,2*J-1] -= (C₁₁₁*(N[i]*B₁[j]+B₁[i]*N[j]) + C₁₂₁*(N[i]*B₂[j]+B₂[i]*N[j]))*ξ.w
                k[2*I-1,2*J]   -= (C₁₂₁*N[i]*B₁[j] + C₁₁₂*B₁[i]*N[j] + C₂₂₁*N[i]*B₂[j] + C₁₂₂*B₂[i]*N[j])*ξ.w
                k[2*I,2*J-1]   -= (C₁₁₂*N[i]*B₁[j] + C₁₂₁*B₁[i]*N[j] + C₁₂₂*N[i]*B₂[j] + C₂₂₁*B₂[i]*N[j])*ξ.w
                k[2*I,2*J]     -= (C₁₂₂*(N[i]*B₁[j]+B₁[i]*N[j]) + C₂₂₂*(N[i]*B₂[j]+B₂[i]*N[j]))*ξ.w
            end
            f[2*I-1] -= ((C₁₁₁*B₁[i]+C₁₂₁*B₂[i])*g₁ + (C₁₁₂*B₁[i]+C₁₂₂*B₂[i])*g₂)*ξ.w
            f[2*I]   -= ((C₁₂₁*B₁[i]+C₂₂₁*B₂[i])*g₁ + (C₁₂₂*B₁[i]+C₂₂₂*B₂[i])*g₂)*ξ.w
        end
    end
end

function (op::Operator{:∫σᵢⱼnⱼgᵢvᵢgᵢds})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N,B₁,B₂ = get∇𝝭(ap,ξ)
        E = op.E
        ν = op.ν
        α = op.α
        Cᵢᵢᵢᵢ = E/(1-ν^2)
        Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
        n₁₁ = ξ.n₁₁
        n₁₂ = ξ.n₁₂
        n₂₂ = ξ.n₂₂
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        C₁₁₁ = Cᵢᵢᵢᵢ*n₁*n₁₁+Cᵢᵢⱼⱼ*n₂*n₁₂
        C₁₁₂ = Cᵢᵢᵢᵢ*n₁*n₁₂+Cᵢᵢⱼⱼ*n₂*n₂₂
        C₂₂₁ = Cᵢᵢⱼⱼ*n₁*n₁₁+Cᵢᵢᵢᵢ*n₂*n₁₂
        C₂₂₂ = Cᵢᵢⱼⱼ*n₁*n₁₂+Cᵢᵢᵢᵢ*n₂*n₂₂
        C₁₂₁ = Cᵢⱼᵢⱼ*(n₁*n₁₂+n₂*n₁₁)
        C₁₂₂ = Cᵢⱼᵢⱼ*(n₂*n₁₂+n₁*n₂₂)
        for i in 1:length(𝓒)
            I = 𝓒[i].id
            for j in 1:length(𝓒)
                J = 𝓒[j].id
                k[2*I-1,2*J-1] -= (C₁₁₁*(N[i]*B₁[j]+B₁[i]*N[j]) + C₁₂₁*(N[i]*B₂[j]+B₂[i]*N[j]) - α*N[i]*n₁₁*N[j])*ξ.w
                k[2*I-1,2*J]   -= (C₁₂₁*N[i]*B₁[j] + C₁₁₂*B₁[i]*N[j] + C₂₂₁*N[i]*B₂[j] + C₁₂₂*B₂[i]*N[j] -α*N[i]*n₁₂*N[j])*ξ.w
                k[2*I,2*J-1]   -= (C₁₁₂*N[i]*B₁[j] + C₁₂₁*B₁[i]*N[j] + C₁₂₂*N[i]*B₂[j] + C₂₂₁*B₂[i]*N[j] -α*N[i]*n₁₂*N[j])*ξ.w
                k[2*I,2*J]     -= (C₁₂₂*(N[i]*B₁[j]+B₁[i]*N[j]) + C₂₂₂*(N[i]*B₂[j]+B₂[i]*N[j]) -α*N[i]*n₂₂*N[j])*ξ.w
            end
            f[2*I-1] -= ((C₁₁₁*B₁[i]+C₁₂₁*B₂[i]-α*N[i]*n₁₁)*g₁ + (C₁₁₂*B₁[i]+C₁₂₂*B₂[i]-α*N[i]*n₁₂)*g₂)*ξ.w
            f[2*I]   -= ((C₁₂₁*B₁[i]+C₂₂₁*B₂[i]-α*N[i]*n₁₂)*g₁ + (C₁₂₂*B₁[i]+C₂₂₂*B₂[i]-α*N[i]*n₂₂)*g₂)*ξ.w
        end
    end
end

function (op::Operator{:∫σ̄ᵢⱼnⱼgᵢds})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N,B₁,B₂ = get∇̄𝝭(ap,ξ)
        E = op.E
        ν = op.ν
        Cᵢᵢᵢᵢ = E/(1-ν^2)
        Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
        n₁₁ = ξ.n₁₁
        n₁₂ = ξ.n₁₂
        n₂₂ = ξ.n₂₂
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        C₁₁₁ = Cᵢᵢᵢᵢ*n₁*n₁₁+Cᵢᵢⱼⱼ*n₂*n₁₂
        C₁₁₂ = Cᵢᵢᵢᵢ*n₁*n₁₂+Cᵢᵢⱼⱼ*n₂*n₂₂
        C₂₂₁ = Cᵢᵢⱼⱼ*n₁*n₁₁+Cᵢᵢᵢᵢ*n₂*n₁₂
        C₂₂₂ = Cᵢᵢⱼⱼ*n₁*n₁₂+Cᵢᵢᵢᵢ*n₂*n₂₂
        C₁₂₁ = Cᵢⱼᵢⱼ*(n₁*n₁₂+n₂*n₁₁)
        C₁₂₂ = Cᵢⱼᵢⱼ*(n₂*n₁₂+n₁*n₂₂)
        for i in 1:length(𝓒)
            I = 𝓒[i].id
            for j in 1:length(𝓒)
                J = 𝓒[j].id
                k[2*I-1,2*J-1] += (C₁₁₁*B₁[i]*N[j] + C₁₂₁*B₂[i]*N[j])*ξ.w
                k[2*I-1,2*J]   += (C₁₁₂*B₁[i]*N[j] + C₁₂₂*B₂[i]*N[j])*ξ.w
                k[2*I,2*J-1]   += (C₁₂₁*B₁[i]*N[j] + C₂₂₁*B₂[i]*N[j])*ξ.w
                k[2*I,2*J]     += (C₁₂₂*B₁[i]*N[j] + C₂₂₂*B₂[i]*N[j])*ξ.w
            end
            f[2*I-1] += ((C₁₁₁*B₁[i] + C₁₂₁*B₂[i])*g₁ + (C₁₁₂*B₁[i] + C₁₂₂*B₂[i])*g₂)*ξ.w
            f[2*I]   += ((C₁₂₁*B₁[i] + C₂₂₁*B₂[i])*g₁ + (C₁₂₂*B₁[i] + C₂₂₂*B₂[i])*g₂)*ξ.w
        end
    end
end

function (op::Operator{:∫σ̃̄ᵢⱼnⱼgᵢds})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N,B₁,B₂ = get∇𝝭(ap,ξ)
        _,B̄₁,B̄₂ = get∇̄𝝭(ap,ξ)
        E = op.E
        ν = op.ν
        Cᵢᵢᵢᵢ = E/(1-ν^2)
        Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
        Cᵢⱼᵢⱼ = E/2/(1+ν)
        n₁₁ = ξ.n₁₁
        n₁₂ = ξ.n₁₂
        n₂₂ = ξ.n₂₂
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        C₁₁₁ = Cᵢᵢᵢᵢ*n₁*n₁₁+Cᵢᵢⱼⱼ*n₂*n₁₂
        C₁₁₂ = Cᵢᵢᵢᵢ*n₁*n₁₂+Cᵢᵢⱼⱼ*n₂*n₂₂
        C₂₂₁ = Cᵢᵢⱼⱼ*n₁*n₁₁+Cᵢᵢᵢᵢ*n₂*n₁₂
        C₂₂₂ = Cᵢᵢⱼⱼ*n₁*n₁₂+Cᵢᵢᵢᵢ*n₂*n₂₂
        C₁₂₁ = Cᵢⱼᵢⱼ*(n₁*n₁₂+n₂*n₁₁)
        C₁₂₂ = Cᵢⱼᵢⱼ*(n₂*n₁₂+n₁*n₂₂)
        for i in 1:length(𝓒)
            I = 𝓒[i].id
            ΔB₁ = B₁[i]-B̄₁[i]
            ΔB₂ = B₂[i]-B̄₂[i]
            for j in 1:length(𝓒)
                J = 𝓒[j].id
                k[2*I-1,2*J-1] -= (C₁₁₁*(N[i]*B₁[j]+ΔB₁*N[j]) + C₁₂₁*(N[i]*B₂[j]+ΔB₂*N[j]))*ξ.w
                k[2*I-1,2*J]   -= (C₁₂₁*N[i]*B₁[j] + C₁₁₂*ΔB₁*N[j] + C₂₂₁*N[i]*B₂[j] + C₁₂₂*ΔB₂*N[j])*ξ.w
                k[2*I,2*J-1]   -= (C₁₁₂*N[i]*B₁[j] + C₁₂₁*ΔB₁*N[j] + C₁₂₂*N[i]*B₂[j] + C₂₂₁*ΔB₂*N[j])*ξ.w
                k[2*I,2*J]     -= (C₁₂₂*(N[i]*B₁[j]+ΔB₁*N[j]) + C₂₂₂*(N[i]*B₂[j]+ΔB₂*N[j]))*ξ.w
            end
            f[2*I-1] -= ((C₁₁₁*ΔB₁+C₁₂₁*ΔB₂)*g₁ + (C₁₁₂*ΔB₁+C₁₂₂*ΔB₂)*g₂)*ξ.w
            f[2*I]   -= ((C₁₂₁*ΔB₁+C₂₂₁*ΔB₂)*g₁ + (C₁₂₂*ΔB₁+C₂₂₂*ΔB₂)*g₂)*ξ.w
        end
    end
end

function (op::Operator{:∫vᵢgᵢds})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = get𝑤(ap,ξ)
        N = get𝝭(ap,ξ)
        n₁₁ = ξ.n₁₁
        n₂₂ = ξ.n₂₂
        n₁₂ = ξ.n₁₂
        g₁ = ξ.g₁
        g₂ = ξ.g₂
        for i in 1:length(𝓒)
            I = 𝓒[i].id
            for j in 1:length(𝓒)
                J = 𝓒[j].id
                k[2*I-1,2*J-1] += op.α*N[i]*n₁₁*N[j]*𝑤
                k[2*I,2*J-1]   += op.α*N[i]*n₁₂*N[j]*𝑤
                k[2*I-1,2*J]   += op.α*N[i]*n₁₂*N[j]*𝑤
                k[2*I,2*J]     += op.α*N[i]*n₂₂*N[j]*𝑤
            end
            f[2*I-1] += op.α*N[i]*(n₁₁*ξ.g₁+n₁₂*ξ.g₂)*𝑤
            f[2*I]   += op.α*N[i]*(n₁₂*ξ.g₁+n₂₂*ξ.g₂)*𝑤
        end
    end
end

## error estimates
function (op::Operator{:L₂})(ap::T) where T<:AbstractElement
    Δu²= 0
    ū² = 0
    for ξ in ap.𝓖
        𝑤 = get𝑤(ap,ξ)
        N = get𝝭(ap,ξ)
        ūᵢ = ξ.u
        uᵢ = 0
        for i in 1:length(ap.𝓒)
            xᵢ = ap.𝓒[i]
            I = xᵢ.id
            uᵢ += N[i]*xᵢ.d
        end
        Δu² += (uᵢ - ūᵢ)^2*𝑤
        ū²  += ūᵢ^2*𝑤
    end
    return Δu², ū²
end

function (op::Operator{:L₂})(aps::Vector{T}) where T<:AbstractElement
    L₂Norm_Δu²= 0
    L₂Norm_ū² = 0
    for ap in aps
        Δu², ū² = op(ap)
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end

function (op::Operator{:H₁})(ap::T) where T<:AbstractElement
    Δ∇u²= 0
    ∇ū² = 0
    Δu²= 0
    ū² = 0
    for ξ in ap.𝓖
        𝑤 = get𝑤(ap,ξ)
        N,B₁,B₂,B₃ = get∇𝝭(ap,ξ)
        ūᵢ = ξ.u
        ∂ūᵢ∂x = ξ.∂u∂x
        ∂ūᵢ∂y = ξ.∂u∂y
        ∂ūᵢ∂z = ξ.∂u∂z
        uᵢ = 0.
        ∂uᵢ∂x = 0.
        ∂uᵢ∂y = 0.
        ∂uᵢ∂z = 0.
        for i in 1:length(ap.𝓒)
            xᵢ = ap.𝓒[i]
            I = xᵢ.id
            uᵢ += N[i]*xᵢ.d
            ∂uᵢ∂x += B₁[i]*xᵢ.d
            ∂uᵢ∂y += B₂[i]*xᵢ.d
            ∂uᵢ∂z += B₃[i]*xᵢ.d
        end
        Δ∇u² += ((∂uᵢ∂x - ∂ūᵢ∂x)^2 + (∂uᵢ∂y - ∂ūᵢ∂y)^2 + (∂uᵢ∂z - ∂ūᵢ∂z)^2)*𝑤
        ∇ū² += (∂ūᵢ∂x^2 + ∂ūᵢ∂y^2 + ∂ūᵢ∂z^2)*𝑤
        Δu² += (uᵢ - ūᵢ)^2*𝑤
        ū² += ūᵢ^2*𝑤
    end
    return Δ∇u², ∇ū², Δu², ū²
end

function (op::Operator{:H₁})(aps::Vector{T}) where T<:AbstractElement
    H₁Norm_Δu²= 0
    H₁Norm_ū² = 0
    L₂Norm_Δu²= 0
    L₂Norm_ū² = 0
    for ap in aps
        Δ∇u², ∇ū², Δu², ū² = op(ap)
        H₁Norm_Δu² += Δu² + Δ∇u²
        H₁Norm_ū²  += ū² + ∇ū²
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (H₁Norm_Δu²/H₁Norm_ū²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end

function (op::Operator{:Hₑ_PlaneStress})(ap::T) where T<:AbstractElement
    ΔW²= 0
    W̄² = 0
    Δu²= 0
    ū² = 0
    E = op.E
    ν = op.ν
    Cᵢᵢᵢᵢ = E/(1-ν^2)
    Cᵢᵢⱼⱼ = E*ν/(1-ν^2)
    Cᵢⱼᵢⱼ = E/2/(1+ν)
    for ξ in ap.𝓖
        N,B₁,B₂ = get∇𝝭(ap,ξ)
        𝑤 = get𝑤(ap,ξ)
        ū₁ = ξ.u
        ū₂ = ξ.v
        ∂ū₁∂x = ξ.∂u∂x
        ∂ū₁∂y = ξ.∂u∂y
        ∂ū₂∂x = ξ.∂v∂x
        ∂ū₂∂y = ξ.∂v∂y
        ε̄₁₁ = ∂ū₁∂x
        ε̄₂₂ = ∂ū₂∂y
        ε̄₁₂ = ∂ū₁∂y + ∂ū₂∂x
        σ̄₁₁ = Cᵢᵢᵢᵢ*ε̄₁₁ + Cᵢᵢⱼⱼ*ε̄₂₂
        σ̄₂₂ = Cᵢᵢᵢᵢ*ε̄₂₂ + Cᵢᵢⱼⱼ*ε̄₁₁
        σ̄₁₂ = Cᵢⱼᵢⱼ*ε̄₁₂
        u₁ = 0.
        u₂ = 0.
        ε₁₁ = 0.
        ε₂₂ = 0.
        ε₁₂ = 0.
        for i in 1:length(ap.𝓒)
            xᵢ = ap.𝓒[i]
            I = xᵢ.id
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            ε₁₁ += B₁[i]*xᵢ.d₁
            ε₂₂ += B₂[i]*xᵢ.d₂
            ε₁₂ += B₂[i]*xᵢ.d₁ + B₁[i]*xᵢ.d₂
        end
        σ₁₁ = Cᵢᵢᵢᵢ*ε₁₁ + Cᵢᵢⱼⱼ*ε₂₂
        σ₂₂ = Cᵢᵢᵢᵢ*ε₂₂ + Cᵢᵢⱼⱼ*ε₁₁
        σ₁₂ = Cᵢⱼᵢⱼ*ε₁₂
        ΔW² += 0.5*((σ₁₁-σ̄₁₁)*(ε₁₁-ε̄₁₁) + (σ₂₂-σ̄₂₂)*(ε₂₂-ε̄₂₂) + (σ₁₂-σ̄₁₂)*(ε₁₂-ε̄₁₂))*𝑤
        W̄² += 0.5*(σ₁₁*ε₁₁ + σ₂₂*ε₂₂ + σ₁₂*ε₁₂)*𝑤
        Δu² += ((u₁ - ū₁)^2 + (u₂ - ū₂)^2)*𝑤
        ū² += (ū₁^2 + ū₂^2)*𝑤
    end
    return ΔW², W̄², Δu², ū²
end

function (op::Operator{:Hₑ_PlaneStress})(aps::Vector{T}) where T<:AbstractElement
    HₑNorm_ΔW²= 0.0
    HₑNorm_W̄² = 0.0
    L₂Norm_Δu²= 0.0
    L₂Norm_ū² = 0.0
    for ap in aps
        ΔW², W̄², Δu², ū² = op(ap)
        HₑNorm_ΔW² += ΔW²
        HₑNorm_W̄²  += W̄²
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end

function set∇𝑢!(ap::T) where T<:AbstractElement
    for ξ in ap.𝓖
        N,B₁,B₂,B₃ = get∇𝝭(ap,ξ)
        𝒙 = get𝒙(ap,ξ)
        u = 0.
        ∂u∂x = 0.
        ∂u∂y = 0.
        ∂u∂z = 0.
        for i in 1:length(ap.𝓒)
            x = ap.𝓒[i]
            I = x.id
            u += N[i]*x.d
            ∂u∂x += B₁[i]*x.d
            ∂u∂y += B₂[i]*x.d
            ∂u∂z += B₃[i]*x.d
        end
        ξ.x = 𝒙[1]
        ξ.y = 𝒙[2]
        ξ.z = 𝒙[3]
        ξ.u = u
        ξ.∂u∂x = ∂u∂x
        ξ.∂u∂y = ∂u∂y
        ξ.∂u∂z = ∂u∂z
    end
end

function (op::Operator{:H₃})(aps::Vector{T}) where T<:AbstractElement
    H₃Norm_Δu²= 0
    H₃Norm_ū² = 0
    H₂Norm_Δu²= 0
    H₂Norm_ū² = 0
    H₁Norm_Δu²= 0
    H₁Norm_ū² = 0
    L₂Norm_Δu²= 0
    L₂Norm_ū² = 0
    for ap in aps
        Δ∇³u², ∇³ū²,Δ∇²u², ∇²ū²,Δ∇u², ∇ū², Δu², ū² = op(ap)
        H₃Norm_Δu² += Δu² + Δ∇u² + Δ∇²u² + Δ∇³u²
        H₃Norm_ū²  += ū² + ∇ū² + ∇²ū² + ∇³ū²
        H₂Norm_Δu² += Δu² + Δ∇u² + Δ∇²u²
        H₂Norm_ū²  += ū² + ∇ū² + ∇²ū²
        H₁Norm_Δu² += Δu² + Δ∇u²
        H₁Norm_ū²  += ū² + ∇ū²
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (H₃Norm_Δu²/H₃Norm_ū²)^0.5, (H₂Norm_Δu²/H₂Norm_ū²)^0.5, (H₁Norm_Δu²/H₁Norm_ū²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end
function (op::Operator{:H₃})(ap::T) where T<:AbstractElement
    Δ∇³u²= 0
    ∇³ū² = 0
    Δ∇²u²= 0
    ∇²ū² = 0
    Δ∇u²= 0
    ∇ū² = 0
    Δu²= 0
    ū² = 0
    for ξ in ap.𝓖
        𝑤 = get𝑤(ap,ξ)
        N,B₁,B₂,B₁₁,B₁₂,B₂₂,B₁₁₁,B₁₁₂,B₁₂₂,B₂₂₂ = get∇³𝝭(ap,ξ)
        ūᵢ = ξ.u
        ∂ūᵢ∂x = ξ.∂u∂x
        ∂ūᵢ∂y = ξ.∂u∂y
        ∂²ūᵢ∂x² = ξ.∂²u∂x²
        ∂²ūᵢ∂x∂y = ξ.∂²u∂x∂y
        ∂²ūᵢ∂y² = ξ.∂²u∂y²
        ∂³ūᵢ∂x³ = ξ.∂³u∂x³
        ∂³ūᵢ∂x²∂y = ξ.∂³u∂x²∂y
        ∂³ūᵢ∂x∂y² = ξ.∂³u∂x∂y²
        ∂³ūᵢ∂y³ = ξ.∂³u∂y³
        uᵢ = 0.
        ∂uᵢ∂x = 0.
        ∂uᵢ∂y = 0.
        ∂²uᵢ∂x² = 0.
        ∂²uᵢ∂x∂y = 0.
        ∂²uᵢ∂y² = 0.
        ∂³uᵢ∂x³ = 0.
        ∂³uᵢ∂x²∂y = 0.
        ∂³uᵢ∂x∂y² = 0.
        ∂³uᵢ∂y³ = 0.
        for i in 1:length(ap.𝓒)
            xᵢ = ap.𝓒[i]
            I = xᵢ.id
            uᵢ += N[i]*xᵢ.d
            ∂uᵢ∂x += B₁[i]*xᵢ.d
            ∂uᵢ∂y += B₂[i]*xᵢ.d
            ∂²uᵢ∂x² += B₁₁[i]*xᵢ.d
            ∂²uᵢ∂x∂y += B₁₂[i]*xᵢ.d
            ∂²uᵢ∂y² += B₂₂[i]*xᵢ.d
            ∂³uᵢ∂x³ += B₁₁₁[i]*xᵢ.d
            ∂³uᵢ∂x²∂y += B₁₁₂[i]*xᵢ.d
            ∂³uᵢ∂x∂y² += B₁₂₂[i]*xᵢ.d
            ∂³uᵢ∂y³ += B₂₂₂[i]*xᵢ.d
        end
        Δ∇³u² += ((∂³uᵢ∂x³ - ∂³ūᵢ∂x³)^2 + (∂³uᵢ∂x²∂y - ∂³ūᵢ∂x²∂y)^2 + (∂³uᵢ∂x∂y² - ∂³ūᵢ∂x∂y²)^2 + (∂³uᵢ∂y³ - ∂³ūᵢ∂y³)^2)*𝑤
        ∇³ū² += (∂³ūᵢ∂x³^2 + ∂³ūᵢ∂x²∂y^2  + ∂³ūᵢ∂x∂y²^2+ ∂³ūᵢ∂y³^2)*𝑤
        Δ∇²u² += ((∂²uᵢ∂x² - ∂²ūᵢ∂x²)^2 + (∂²uᵢ∂x∂y - ∂²ūᵢ∂x∂y)^2 + (∂²uᵢ∂y² - ∂²ūᵢ∂y²)^2)*𝑤
        ∇²ū² += (∂²ūᵢ∂x²^2 + ∂²ūᵢ∂x∂y^2 + ∂²ūᵢ∂y²^2)*𝑤
        Δ∇u² += ((∂uᵢ∂x - ∂ūᵢ∂x)^2 + (∂uᵢ∂y - ∂ūᵢ∂y)^2)*𝑤
        ∇ū² += (∂ūᵢ∂x^2 + ∂ūᵢ∂y^2)*𝑤
        Δu² += (uᵢ - ūᵢ)^2*𝑤
        ū² += ūᵢ^2*𝑤
    end
    return Δ∇³u², ∇³ū², Δ∇²u², ∇²ū², Δ∇u², ∇ū², Δu², ū²
end

function set∇𝑢!(aps::Vector{T}) where T<:AbstractElement
    for ap in aps
        set∇𝑢!(ap)
    end
end
