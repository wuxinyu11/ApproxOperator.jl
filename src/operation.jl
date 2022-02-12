
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

@inline function (op::Operator)(aps::Vector{T},gps::Vector{S},k::AbstractMatrix{Float64},f::Vector{Float64}) where {T<:AbstractElement,S<:AbstractElement}
    for i in 1:length(aps)
        op(aps[i],gps[i],k,f)
    end
end

@inline function (op::Operator)(aps::Vector{T},k::AbstractMatrix{Float64},f::Vector{Float64}) where T<:AbstractElement
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
        v = f(𝒙...)
        setproperty!(ξ,s,v)
    end
end

function prescribe!(aps::Vector{T},s::Symbol,f::Function) where T<:AbstractElement
    𝓖 = aps[1].𝓖
    data = 𝓖[1].data
    if ~haskey(data,s)
        push!(data,s=>similar(data[:w]))
    end
    for ap in aps
        prescribe!(ap,s,f)
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
        ~,B₁,B₂,B₃ = get∇𝝭(ap,ξ)
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
                k[I,J] -= ((B₁[i]*n₁+B₂[i]*n₂+B₃[i]*n₃)*N[j] + N[i]*(B₁[j]*n₁+B₂[j]*n₂+B₃[j]*n₃))*ξ.w
            end
            f[I] -= (B₁[i]*n₁+B₂[i]*n₂+B₃[i]*n₃)*ξ.g*ξ.w
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
        end
    end
end

function (op::Operator{:∫∫vᵢbᵢdxdy})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
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
        for i in 1:length(𝓒)
            I = 𝓒[i].id
            for j in 1:length(𝓒)
                J = 𝓒[j].id
                k[2*I-1,2*J-1] -= ((Cᵢᵢᵢᵢ*n₁*n₁₁+Cᵢᵢⱼⱼ*n₂*n₁₂)*(N[i]*B₁[j]+B₁[i]*N[j]) + Cᵢⱼᵢⱼ*(n₁*n₁₂+n₂*n₁₁)*(N[i]*B₂[j]+B₂[i]*N[j]))*ξ.w
                k[2*I-1,2*J]   -= (Cᵢⱼᵢⱼ*(n₁*n₁₂+n₂*n₁₁)*N[i]*B₁[j] + (Cᵢᵢᵢᵢ*n₁*n₁₂+Cᵢᵢⱼⱼ*n₂*n₂₂)*B₁[i]*N[j] + (Cᵢᵢⱼⱼ*n₁*n₁₁+Cᵢᵢᵢᵢ*n₂*n₁₂)*N[i]*B₂[j] + Cᵢⱼᵢⱼ*(n₂*n₁₂+n₁*n₂₂)*B₂[i]*N[j])*ξ.w
                k[2*I,2*J-1]   -= ((Cᵢᵢᵢᵢ*n₁*n₁₂+Cᵢᵢⱼⱼ*n₂*n₂₂)*N[i]*B₁[j] + Cᵢⱼᵢⱼ*(n₁*n₁₂+n₂*n₁₁)*B₁[i]*N[j] + Cᵢⱼᵢⱼ*(n₂*n₁₂+n₁*n₂₂)*N[i]*B₂[j] + (Cᵢᵢⱼⱼ*n₁*n₁₁+Cᵢᵢᵢᵢ*n₂*n₁₂)*B₂[i]*N[j])*ξ.w
                k[2*I,2*J]     -= (Cᵢⱼᵢⱼ*(n₂*n₁₂+n₁*n₂₂)*(N[i]*B₁[j]+B₁[i]*N[j]) + (Cᵢᵢⱼⱼ*n₁*n₁₂+Cᵢᵢᵢᵢ*n₂*n₂₂)*(N[i]*B₂[j]+B₂[i]*N[j]))*ξ.w
            end
            f[2*I-1] -= (((Cᵢᵢᵢᵢ*n₁*n₁₁+Cᵢᵢⱼⱼ*n₂*n₁₂)*B₁[i] + Cᵢⱼᵢⱼ*(n₁*n₁₂+n₂*n₁₁)*B₂[i])*g₁ + ((Cᵢᵢᵢᵢ*n₁*n₁₂+Cᵢᵢⱼⱼ*n₂*n₂₂)*B₁[i] + Cᵢⱼᵢⱼ*(n₂*n₁₂+n₁*n₂₂)*B₂[i])*g₂)*ξ.w
            f[2*I]   -= ((Cᵢⱼᵢⱼ*(n₁*n₁₂+n₂*n₁₁)*B₁[i] + (Cᵢᵢⱼⱼ*n₁*n₁₁+Cᵢᵢᵢᵢ*n₂*n₁₂)*B₂[i])*g₁ + (Cᵢⱼᵢⱼ*(n₂*n₁₂+n₁*n₂₂)*B₁[i] + (Cᵢᵢⱼⱼ*n₁*n₁₂+Cᵢᵢᵢᵢ*n₂*n₂₂)*B₂[i])*g₂)*ξ.w
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
    HₑNorm_ΔW²= 0
    HₑNorm_W̄² = 0
    L₂Norm_Δu²= 0
    L₂Norm_ū² = 0
    for ap in aps
        ΔW², W̄², Δu², ū² = op(ap)
        HₑNorm_ΔW² += ΔW²
        HₑNorm_W̄²  += W̄²
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (HₑNorm_ΔW²/HₑNorm_W̄²)^0.5, (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end
