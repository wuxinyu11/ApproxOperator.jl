## Counstruction
struct Operator{T,D}
    type::Val{T}
    data::Dict{Symbol,D}
end
Operator(type::Symbol,data::Dict{Symbol,D}) where D = Operator(Val(type),data)
Operator(type::Symbol) = Operator(Val(type),Dict{Symbol,Float64}())

## General Functions
push!(op::Operator{T,D},d::Pair{Symbol,D}...) where {T,D} = push!(op.data,d...)
@inline getproperty(op::Operator,f::Symbol) = hasfield(Operator,f) ? getfield(op,f) : getfield(op,:data)[f]
@inline function (op::Operator)(aps::Vector{T},k::AbstractMatrix{Float64},f::Vector{Float64}) where T<:Approximator
    for ap in aps
        op(ap,k,f)
    end
end
@inline function (op::Operator)(aps::Vector{T},f::Vector{Float64}) where T<:Approximator
    for ap in aps
        op(ap,f)
    end
end

function (op::Operator{:update})(ap::Approximator,s::Symbol,f::Function)
    𝓖 = ap.𝓖
    for ξ in 𝓖
        x = getx(ap,ξ)
        v = f(x...)
        ξ.s = v
    end
end

## Potential Problem
function (op::Operator{:∫∇v∇udΩ})(ap::Approximator,k::AbstractMatrix{Float64},f::AbstractVector{Float64})
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N,B₁,B₂,B₃ = get∇𝝭(ap,ξ)
        w = getw(ap,ξ)
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

function (op::Operator{:∫vtdΓ})(ap::Approximator,f::AbstractVector{Float64})
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        w = getw(ap,ξ)
        N = get𝝭(ap,ξ)
        for i in 1:length(𝓒)
            I = 𝓒[i].id
            f[I] = f[I] + N[i]*ξ.t*w
        end
    end
end

function (op::Operator{:g})(ap::Poi1,k::AbstractMatrix{Float64},f::AbstractVector{Float64};dof::Symbol=:d)
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
