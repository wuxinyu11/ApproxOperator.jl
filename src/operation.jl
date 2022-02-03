
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

@inline function (op::Operator)(aps::Vector{T},gps::Vector{S},k::AbstractMatrix{Float64},f::Vector{Float64}) where {T<:Approximator,S<:Approximator}
    for i in 1:length(aps)
        op(aps[i],gps[i],k,f)
    end
end

@inline function (op::Operator)(aps::Vector{T},k::AbstractMatrix{Float64},f::Vector{Float64}) where T<:Approximator
    for ap in aps
        op(ap,k,f)
    end
end
@inline function (op::Operator)(aps::Vector{T},k::AbstractMatrix{Float64}) where T<:Approximator
    for ap in aps
        op(ap,k)
    end
end
@inline function (op::Operator)(aps::Vector{T},f::Vector{Float64}) where T<:Approximator
    for ap in aps
        op(ap,f)
    end
end

@inline function (op::Operator)(aps::Vector{T},s::Symbol) where T<:Approximator
    for ap in aps
        op(ap,s)
    end
end
@inline function (op::Operator)(aps::Vector{T}) where T<:Approximator
    for ap in aps
        op(ap)
    end
end

function prescribe!(ap::T,s::Symbol,f::Function) where T<:Approximator
    𝓖 = ap.𝓖
    for ξ in 𝓖
        x = getx(ap,ξ)
        v = f(x...)
        setproperty!(ξ,s,v)
    end
end

function prescribe!(aps::Vector{T},s::Symbol,f::Function) where T<:Approximator
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
function (op::Operator{:∫∇v∇uvbdΩ})(ap::Approximator,k::AbstractMatrix{Float64},f::AbstractVector{Float64})
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

function (op::Operator{:∫∇v∇udΩ})(ap::Approximator,k::AbstractMatrix{Float64})
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

function (op::Operator{:∫vbdΩ})(ap::Approximator,f::AbstractVector{Float64})
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = get𝝭(ap,ξ)
        𝑤 = get𝑤(ap,ξ)
        for i in 1:length(𝓒)
            I = 𝓒[i].id
            f[I] += N[i]*ξ.b*𝑤
        end
    end
end

function (op::Operator{:∫vtdΓ})(ap::Approximator,f::AbstractVector{Float64})
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

function (op::Operator{:∫vgdΓ})(ap::Approximator,k::AbstractMatrix{Float64},f::AbstractVector{Float64})
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

function (op::Operator{:∫λgdΓ})(ap1::Approximator,ap2::Approximator,g::AbstractMatrix{Float64},q::AbstractVector{Float64})
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

function (op::Operator{:∫∇𝑛vgdΓ})(ap::Approximator,k::AbstractMatrix{Float64},f::AbstractVector{Float64})
    for ξ in ap.𝓖
        𝑤 = get𝑤(ap,ξ)
        N,B = get∇𝑛𝝭(ap,ξ)
        for i in 1:length(ap.𝓒)
            I = ap.𝓒[i].id
            for j in 1:length(ap.𝓒)
                J = ap.𝓒[j].id
                k[I,J] += (-B[i]*N[j] - N[i]*B[j] + op.α*N[i]*N[j])*𝑤
            end
            f[I] += (op.α*N[i]*ξ.g̃ - B[i]*ξ.g)*𝑤
        end
    end
end

function (op::Operator{:g})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64};dof::Symbol=:d) where T<:AbstractPoi
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

## error estimates
function (op::Operator{:L₂})(ap::T) where T<:Approximator
    Δu²= 0
    ū² = 0
    for ξ in ap.𝓖
        w = getw(ap,ξ)
        N = get𝝭(ap,ξ)
        x = getx(ap,ξ)
        ūᵢ = ξ.u
        uᵢ = 0
        for i in 1:length(ap.𝓒)
            xᵢ = ap.𝓒[i]
            I = xᵢ.id
            uᵢ += N[i]*xᵢ.d
        end
        Δu² += (uᵢ - ūᵢ)^2*w
        ū²  += ūᵢ^2*w
    end
    return Δu², ū²
end

function (op::Operator{:L₂})(aps::Vector{T}) where T<:Approximator
    L₂Norm_Δu²= 0
    L₂Norm_ū² = 0
    for ap in aps
        Δu², ū² = op(ap)
        L₂Norm_Δu² += Δu²
        L₂Norm_ū²  += ū²
    end
    return (L₂Norm_Δu²/L₂Norm_ū²)^0.5
end

function (op::Operator{:H₁})(ap::T) where T<:Approximator
    Δ∇u²= 0
    ∇ū² = 0
    Δu²= 0
    ū² = 0
    for ξ in ap.𝓖
        w = getw(ap,ξ)
        N,B₁,B₂,B₃ = get∇𝝭(ap,ξ)
        x = getx(ap,ξ)
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
        Δ∇u² += ((∂uᵢ∂x - ∂ūᵢ∂x)^2 + (∂uᵢ∂y - ∂ūᵢ∂y)^2 + (∂uᵢ∂z - ∂ūᵢ∂z)^2)*w
        ∇ū² += (∂ūᵢ∂x^2 + ∂ūᵢ∂y^2 + ∂ūᵢ∂z^2)*w
        Δu² += (uᵢ - ūᵢ)^2*w
        ū² += ūᵢ^2*w
    end
    return Δ∇u², ∇ū², Δu², ū²
end

function (op::Operator{:H₁})(aps::Vector{T}) where T<:Approximator
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
