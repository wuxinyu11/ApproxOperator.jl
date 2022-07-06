
# function checkIC(a::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Tri3},b::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Tri3}) where {𝒑,𝑠,𝜙}
#     x₁ = a.𝓒[1].x;y₁ = a.𝓒[1].y
#     x₂ = a.𝓒[2].x;y₂ = a.𝓒[2].y
#     x₃ = a.𝓒[3].x;y₃ = a.𝓒[3].y
#     n₁₁ = y₃-y₂;n₂₁ = y₁-y₃;n₃₁ = y₂-y₁
#     n₁₂ = x₂-x₃;n₂₂ = x₃-x₁;n₃₂ = x₁-x₂
#     nᵈ = length(a.𝓒)
#     nᵖ = get𝑛𝒒(a)
#     # nᵖ = get𝑛𝒑(a)
#     f₁ = zeros(nᵈ,nᵖ)
#     f₂ = zeros(nᵈ,nᵖ)
#     for ξ in a.𝓖
#         N,B₁,B₂ = get∇𝝭(a,ξ)
#         𝒙 = get𝒙(a,ξ)
#         𝒒 = get𝒒(a,ξ)
#         p = get𝒑(a,𝒙)
#         𝑤 = get𝑤(a,ξ)
#         for i in 1:nᵈ
#             for j in 1:nᵖ
#                 f₁[i,j] += B₁[i]*p[j]*𝑤
#                 f₂[i,j] += B₂[i]*p[j]*𝑤
#             end
#         end
#     end
#     for ξ in b.𝓖
#         N = get𝝭(b,ξ)
#         𝒒,∂𝒒∂ξ,∂𝒒∂η = get∇𝒒(b,ξ)
#         𝒙 = get𝒙(b,ξ)
#         p,∂p∂x,∂p∂y = get∇𝒑(b,𝒙)
#         wᵇ = ξ.wᵇ
#         w = ξ.w
#         𝑤 = get𝑤(b,ξ)
#         nᵇ₁ = 0.0;nᵇ₂ = 0.0
#         ξ.ξ == 0.0 ? (nᵇ₁ += n₁₁;nᵇ₂ += n₁₂) : nothing
#         ξ.η == 0.0 ? (nᵇ₁ += n₂₁;nᵇ₂ += n₂₂) : nothing
#         ξ.ξ+ξ.η ≈ 1.0 ? (nᵇ₁ += n₃₁;nᵇ₂ += n₃₂) : nothing
#         for j in 1:nᵖ
#             W₁ = p[j]*nᵇ₁*wᵇ - ∂p∂x[j]*𝑤
#             W₂ = p[j]*nᵇ₂*wᵇ - ∂p∂y[j]*𝑤
#             for i in 1:nᵈ
#                 f₁[i,j] -= N[i]*W₁
#                 f₂[i,j] -= N[i]*W₂
#             end
#         end
#     end
#     return f₁,f₂
# end

# function checkCC(a::ReproducingKernel)
#     nᵖ = get𝑛𝒑(a)
#     nⁱ = length(a.𝓖)
#     f = zeros(nⁱ,nᵖ)
#     f₁ = zeros(nⁱ,nᵖ)
#     f₂ = zeros(nⁱ,nᵖ)
#     for i in 1:length(a.𝓖)
#         ξ = a.𝓖[i]
#         𝒙 = get𝒙(a,ξ)
#         p,∂p∂x,∂p∂y = get∇𝒑(a,𝒙)
#         N,B₁,B₂ = get∇𝝭(a,ξ)
#         for j in 1:nᵖ
#             for k in 1:length(a.𝓒)
#                 𝒙̄ = (a.𝓒[k].x,a.𝓒[k].y,a.𝓒[k].z)
#                 p̄ = get𝒑(a,𝒙̄)
#                 f[i,j]  += N[k]*p̄[j]
#                 f₁[i,j] += B₁[k]*p̄[j]
#                 f₂[i,j] += B₂[k]*p̄[j]
#             end
#             f[i,j]  -= p[j]
#             f₁[i,j] -= ∂p∂x[j]
#             f₂[i,j] -= ∂p∂y[j]
#         end
#     end
#     return f,f₁,f₂
# end

function checkConsistency(as::Vector{T}) where T<:ReproducingKernel
    nᵖ = get𝑛𝒑(as[1])
    f = zeros(nᵖ)
    𝒑ʰ = zeros(nᵖ)
    set𝝭!(as)
    for a in as
        checkConsistency(a,f,𝒑ʰ)
    end
    return f.^0.5
end

const list∇𝝭 = (:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z)
const list∇²𝝭 = (:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂²𝝭∂x²,:∂²𝝭∂x∂y,:∂²𝝭∂y²,:∂𝝭∂z,:∂²𝝭∂x∂z,:∂²𝝭∂y∂z,:∂²𝝭∂z²)
const list∇³𝝭 = (:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂²𝝭∂x²,:∂²𝝭∂x∂y,:∂²𝝭∂y²,:∂³𝝭∂x³,:∂³𝝭∂x²∂y,:∂³𝝭∂x∂y²,:∂³𝝭∂y³)
for (𝝭,𝒑,list) in ((:check𝝭,:get𝒑,:list𝝭),
                   (:check∇𝝭,:get∇𝒑,:list∇𝝭),
                   (:check∇²𝝭,:get∇²𝒑,:list∇²𝝭),
                   (:check∇³𝝭,:get∇³𝒑,:list∇³𝝭))
    @eval begin
        function $𝝭(a::T,f::Matrix{Float64},𝒑::Matrix{Float64},𝒑ʰ::Matrix{Float64}) where T<:AbstractElement
            n = get𝑛𝒑(a)
            for ξ in a.𝓖
                𝒙 = get𝒙(a,ξ)
                𝑤 = ξ.𝑤
                𝒑s = $𝒑(a,𝒙)
                for i in 1:n
                    for (j,𝒑_) in enumerate(𝒑s)
                        𝒑[i,j] = 𝒑_[i]
                    end
                end
                fill!(𝒑ʰ,0.0)
                for (k,𝒙ᵢ) in enumerate(a.𝓒)
                    𝒑ᵢ = get𝒑(a,(𝒙ᵢ.x,𝒙ᵢ.y,𝒙ᵢ.z))
                    for i in 1:n
                        for (j,s) in enumerate($list)
                            𝒑ʰ[i,j] += ξ[s][k]*𝒑ᵢ[i]
                        end
                    end
                end
                f .+= (𝒑 .- 𝒑ʰ).^2 .* 𝑤
            end
        end

        function $𝝭(as::Vector{T}) where T<:ReproducingKernel
            nᵖ = get𝑛𝒑(as[1])
            n = length($list)
            f = zeros(nᵖ,n)
            𝒑 = zeros(nᵖ,n)
            𝒑ʰ = zeros(nᵖ,n)
            for a in as
                $𝝭(a,f,𝒑,𝒑ʰ)
            end
            return f.^0.5
        end
    end
end
