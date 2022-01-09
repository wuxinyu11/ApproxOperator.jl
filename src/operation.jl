
struct Operator{T<:Any}
    type::Val{T}
    data::Dict{Symbol,Float64}
end
@inline getproperty(op::Operator,f::Symbol) = hasfield(Operator,f) ? getfield(op,f) : op.data[f]
# @inline (op::Operator)(ap::T,k::AbstractMatrix{Float64},f::Vector{Float64}) where T<:Approximator = op(ap,k,f,Val(op.type))
# @inline (op::Operator)(ap::T,f::Vector{Float64}) where T<:Approximator = op(ap,f,op.type)
# @inline function (op::Operator)(aps::Vector{T},k::AbstractMatrix{Float64},f::Vector{Float64}) where T<:Approximator
#     for ap in aps
#         op(ap,k,f,op.type)
#     end
# end
# @inline function (op::Operator)(aps::Vector{T},f::Vector{Float64}) where T<:Approximator
#     for ap in aps
#         op(ap,f,op.type)
#     end
# end


## Potential Problem
function (op::Operator{:∇v∇u})(ap::Approximator,k::AbstractMatrix{Float64},f::AbstractVector{Float64})
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

function (op::Operator{:vt})(ap::Approximator,f::AbstractVector{Float64})
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
