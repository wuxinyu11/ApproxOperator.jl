
struct Operator{T<:Any}
    type::Val{T}
    data::Dict{Symbol,Float64}
end
Operator(type::Symbol) = Operator(Val(type),Dict{Symbol,Float64}())

@inline getproperty(op::Operator,f::Symbol) = hasfield(Operator,f) ? getfield(op,f) : op.data[f]
@inline function (op::Operator{S})(aps::Vector{T},k::AbstractMatrix{Float64},f::Vector{Float64}) where {S<:Any,T<:Approximator}
    for ap in aps
        op(ap,k,f)
    end
end
@inline function (op::Operator{S})(aps::Vector{T},f::Vector{Float64}) where {S<:Any,T<:Approximator}
    for ap in aps
        op(ap,f)
    end
end

function (op::Operator{:∫∇v∇udΩ})(ap::Approximator,k::AbstractMatrix{Float64},f::AbstractVector{Float64})
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N,B₁,B₂,B₃ = ap.∇𝝭(ξ)
        w = ap.J(ξ)*ξ.w
        x = ap.coordinates(ξ)
        for i in 1:length(𝓒)
            I = 𝓒[i].id
            for j in 1:length(𝓒)
                J = 𝓒[j].id
                k[I,J] += op.k*(B₁[i]*B₁[j] + B₂[i]*B₂[j] + B₃[i]*B₃[j])*w
            end
            f[I] += N[i]*ξ.b*w
        end
    end
end

## Potential Problem
function (op::Operator{:∫∇v∇udΩ})(ap::Approximator,k::AbstractMatrix{Float64},f::AbstractVector{Float64})
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = ap.𝝭(ξ)
        B₁,B₂,B₃ = ap.∇𝝭(ξ)
        w = ap.J(ξ)*ξ.w
        x = ap.coordinates(ξ)
        for i in 1:length(𝓒)
            I = 𝓒[i].I
            for j in 1:length(𝓒)
                J = 𝓒[j].I
                k[I,J] += op.k*(B₁[i]*B₁[j] + B₂[i]*B₂[j] + B₃[i]*B₃[j])*w
            end
            f[I] += N[i]*ξ.b*w
        end
    end
end

function (op::Operator{:∫vtdΓ})(ap::Approximator,f::AbstractVector{Float64})
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        w = ap.J(ξ)*ξ.w
        N = ap.𝝭(ξ)
        x = ap.coordinates(ξ)
        for i in 1:length(𝓒)
            I = 𝓒[i].I
            f[I] = f[I] + N[i]*ξ.t*w
        end
    end
end

function (op::Operator{:g})(x::Node,k::AbstractMatrix{Float64},f::AbstractVector{Float64})
    j = x.I
    for i in 1:length(f)
        f[i] -= k[i,j]*x.d
    end
    k[j,:] .= 0.
    k[:,j] .= 0.
    k[j,j] = 1.
    f[j] = x.g
end
