
const list𝝭 = (:𝝭,)
const list∇𝝭 = (:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z)
const list∇₁𝝭 = (:𝝭,:∂𝝭∂x)
const list∇₂𝝭 = (:𝝭,:∂𝝭∂x,:∂𝝭∂y)
const list∇²𝝭 = (:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂²𝝭∂x²,:∂²𝝭∂x∂y,:∂²𝝭∂y²,:∂𝝭∂z,:∂²𝝭∂x∂z,:∂²𝝭∂y∂z,:∂²𝝭∂z²)
const list∇²₂𝝭 = (:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂²𝝭∂x²,:∂²𝝭∂x∂y,:∂²𝝭∂y²)
const list∇̃²₂𝝭 = (:∂²𝝭∂x²,:∂²𝝭∂x∂y,:∂²𝝭∂y∂x,:∂²𝝭∂y²)
const list∇³𝝭 = (:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂²𝝭∂x²,:∂²𝝭∂x∂y,:∂²𝝭∂y²,:∂³𝝭∂x³,:∂³𝝭∂x²∂y,:∂³𝝭∂x∂y²,:∂³𝝭∂y³)
const list∇∇²𝝭 = (:𝝭,:∂∂²𝝭∂x²∂x,:∂∂²𝝭∂x²∂y,:∂∂²𝝭∂x∂y∂x,:∂∂²𝝭∂x∂y∂y,:∂∂²𝝭∂y²∂x,:∂∂²𝝭∂y²∂y)
for (𝝭,𝒑,list) in ((:check∇𝝭,:get∇𝒑,:list∇𝝭),
                   (:check∇₁𝝭,:get∇₁𝒑,:list∇₁𝝭),
                   (:check∇₂𝝭,:get∇₂𝒑,:list∇₂𝝭),
                   (:check∇²𝝭,:get∇²𝒑,:list∇²𝝭),
                   (:check∇∇²𝝭,:get∇∇²𝒑,:list∇∇²𝝭),
                   (:check∇²₂𝝭,:get∇²₂𝒑,:list∇²₂𝝭),
                   (:check∇̃²₂𝝭,:get∇̃²₂𝒑,:list∇̃²₂𝝭),
                   (:check∇³𝝭,:get∇³𝒑,:list∇³𝝭))
    @eval begin
        function $𝝭(a::T,f::Matrix{Float64},𝒑::Matrix{Float64},𝒑ʰ::Matrix{Float64}) where T<:AbstractElement
            n = get𝑛𝒑(a)
            for ξ in a.𝓖
                𝑤 = ξ.𝑤
                𝒑s = $𝒑(a,(ξ.x,ξ.y,ξ.z))
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

function check𝝭(a::T,f::Vector{Float64},𝒑::Vector{Float64},𝒑ʰ::Vector{Float64}) where T<:AbstractElement
    n = get𝑛𝒑(a)
    for ξ in a.𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        𝒑 = get𝒑(a,(ξ.x,ξ.y,ξ.z))
        fill!(𝒑ʰ,0.0)
        for (k,𝒙ᵢ) in enumerate(a.𝓒)
            𝒑ᵢ = get𝒑(a,(𝒙ᵢ.x,𝒙ᵢ.y,𝒙ᵢ.z))
            for i in 1:n
                𝒑ʰ[i] += N[k]*𝒑ᵢ[i]
            end
        end
        f .+= (𝒑 .- 𝒑ʰ).^2 .* 𝑤
    end
end
    
function check𝝭(as::Vector{T}) where T<:ReproducingKernel
    nᵖ = get𝑛𝒑(as[1])
    f = zeros(nᵖ)
    𝒑 = zeros(nᵖ)
    𝒑ʰ = zeros(nᵖ)
    for a in as
        check𝝭(a,f,𝒑,𝒑ʰ)
    end
    return f.^0.5
end

# function check𝝭ₕ(as::Vector{T}) where T<:ReproducingKernel
#     nᵖ = get𝑛𝒑(as[1])
#     f = zeros(nᵖ)
#     𝒑 = zeros(nᵖ)
#     𝒑ʰ = zeros(nᵖ)
#     for a in as
#         check𝝭ₕ(a,f,𝒑,𝒑ʰ)
#     end
#     return f.^0.5
# end

# function check𝝭(a::ReproducingKernel{:Linear2D},f::Vector{Float64},𝒑::Vector{Float64},𝒑ʰ::Vector{Float64})
#     n = get𝑛𝒑(a)
#     for ξ in a.𝓖
#         𝑤 = ξ.𝑤
#         N = ξ[:𝝭]
#         N₁ = ξ[:𝝭₁]
#         N₂ = ξ[:𝝭₂]
#         𝒑 = get𝒑(a,(ξ.x,ξ.y,ξ.z))
#         fill!(𝒑ʰ,0.0)
#         for (k,𝒙ᵢ) in enumerate(a.𝓒)
#             # 𝒑ᵢ,∂𝒑ᵢ∂x,∂𝒑ᵢ∂y = get∇₂𝒑(a,(𝒙ᵢ.x,𝒙ᵢ.y,𝒙ᵢ.z))
#             # 𝒑ᵢ,∂𝒑ᵢ∂x,∂𝒑ᵢ∂y = get∇₂𝒑(a,(𝒙ᵢ.x-ξ.x,𝒙ᵢ.y-ξ.y,𝒙ᵢ.z-ξ.z))
#             𝒑ᵢ,∂𝒑ᵢ∂x,∂𝒑ᵢ∂y = get∇₂𝒑(a,(ξ.x-𝒙ᵢ.x,ξ.y-𝒙ᵢ.y,ξ.z-𝒙ᵢ.z))
#             for i in 1:n
#                 𝒑ʰ[i] += N[k]*𝒑ᵢ[i] + N₁[k]*∂𝒑ᵢ∂x[i] + N₂[k]*∂𝒑ᵢ∂y[i]
#             end
#         end
#         f .+= (𝒑 .- 𝒑ʰ).^2 .* 𝑤
#     end
# end
function check∇₂𝝭(a::ReproducingKernel{:Linear2D},f::Matrix{Float64},𝒑::Matrix{Float64},𝒑ʰ::Matrix{Float64})
    n = get𝑛𝒑(a)
    for ξ in a.𝓖
        𝑤 = ξ.𝑤
        𝒑s = get∇₂𝒑(a,(0.0,0.0,0.0))
        for i in 1:n
            for (j,𝒑_) in enumerate(𝒑s)
                𝒑[i,j] = 𝒑_[i]
            end
        end
        fill!(𝒑ʰ,0.0)
        for (k,𝒙ᵢ) in enumerate(a.𝓒)
            𝒑ᵢ,∂𝒑ᵢ∂x,∂𝒑ᵢ∂y = get∇₂𝒑(a,(ξ.x-𝒙ᵢ.x,ξ.y-𝒙ᵢ.y,ξ.z-𝒙ᵢ.z))
            for i in 1:n
                for (j,s) in enumerate(((:𝝭,:𝝭₁,:𝝭₂),(:∂𝝭∂x,:∂𝝭₁∂x,:∂𝝭₂∂x),(:∂𝝭∂y,:∂𝝭₁∂y,:∂𝝭₂∂y)))
                    𝒑ʰ[i,j] += ξ[s[1]][k]*𝒑ᵢ[i] + ξ[s[2]][k]*∂𝒑ᵢ∂x[i] + ξ[s[3]][k]*∂𝒑ᵢ∂y[i]
                end
            end
        end
        f .+= (𝒑 .- 𝒑ʰ).^2 .* 𝑤
    end
end