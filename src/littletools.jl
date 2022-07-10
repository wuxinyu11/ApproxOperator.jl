
const list𝝭 = (:𝝭,)
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
