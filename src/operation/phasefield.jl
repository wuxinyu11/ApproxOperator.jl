
function (op::Operator{:∫v²uₓuₓdx})(ap::T,k::AbstractMatrix{Float64},f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    EA = op.EA
    η = op.η
    for ξ in 𝓖
        N = ξ[:𝝭]
        B = ξ[:∂𝝭∂x]
        v = sum(N[i]*xᵢ.v for (i,xᵢ) in enumerate(𝓒))
        # println(v)
        b = ξ.b
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                k[I,J] += (v^2+η)*EA*B[i]*B[j]*𝑤
            end
            f[I] += N[i]*b*𝑤
        end
    end
end
