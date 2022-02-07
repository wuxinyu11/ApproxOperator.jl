
## ReproducingKernel
struct ReproducingKernel{𝝃,𝑝,𝑠,𝜙,T}<:AbstractElement{T}
    𝓒::Vector{Node}
    𝓖::Vector{𝝃}
    𝗠::Dict{Symbol,SymMat}
    𝝭::Dict{Symbol,Vector{Float64}}
end

## shape functions
function get𝝭(ap::ReproducingKernel,ξ::Node)
    𝓒 = ap.𝓒
    𝝭 = ap.𝝭[:∂1]
    𝒙 = get𝒙(ap,ξ)
    𝒑₀ᵀ𝗠⁻¹ = cal𝗠!(ap,𝒙)
    for i in 1:length(𝓒)
        𝒙ᵢ = 𝓒[i]
        Δ𝒙 = 𝒙 - 𝒙ᵢ
        𝒑 = get𝒑(ap,Δ𝒙)
        𝜙 = get𝜙(ap,𝒙ᵢ,Δ𝒙)
        𝝭[i] = 𝒑₀ᵀ𝗠⁻¹*𝒑*𝜙
    end
    return 𝝭
end

function get∂𝝭∂x(ap::ReproducingKernel,ξ::Node)
    𝓒 = ap.𝓒
    𝝭 = ap.𝝭[:∂1]
    ∂𝝭∂x = ap.𝝭[:∂x]
    𝒙 = getx(ap,ξ)
    𝒑₀ᵀ𝗠⁻¹, 𝒑₀ᵀ∂𝗠⁻¹∂x = cal∂𝗠∂x!(ap,𝒙)
    for i in 1:length(𝓒)
        𝒙ᵢ = 𝓒[i]
        Δ𝒙 = 𝒙 - 𝒙ᵢ
        𝒑, ∂𝒑∂x = get∂𝒑∂x(ap,Δx)
        𝜙, ∂𝜙∂x = get∂𝜙∂x(ap,xᵢ,Δx)
        𝝭[i] = 𝒑₀ᵀ𝗠⁻¹*𝒑*𝜙
        ∂𝝭∂x[i] = 𝒑₀ᵀ∂𝗠⁻¹∂x*𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂x*𝜙 + 𝒑₀ᵀ𝗠⁻¹*𝒑*∂𝜙∂x
    end
    return 𝝭, ∂𝝭∂x
end

function get∇𝝭(ap::ReproducingKernel,ξ::Node)
    𝓒 = ap.𝓒
    𝝭 = ap.𝝭[:∂1]
    ∂𝝭∂x = ap.𝝭[:∂x]
    ∂𝝭∂y = ap.𝝭[:∂y]
    ∂𝝭∂z = ap.𝝭[:∂z]
    𝒙 = get𝒙(ap,ξ)
    𝒑₀ᵀ𝗠⁻¹, 𝒑₀ᵀ∂𝗠⁻¹∂x, 𝒑₀ᵀ∂𝗠⁻¹∂y, 𝒑₀ᵀ∂𝗠⁻¹∂z= cal∇𝗠!(ap,𝒙)
    for i in 1:length(𝓒)
        𝒙ᵢ = 𝓒[i]
        Δ𝒙 = 𝒙 - 𝒙ᵢ
        𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂𝒑∂z = get∇𝒑(ap,Δ𝒙)
        𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂𝜙∂z = get∇𝜙(ap,𝒙ᵢ,Δ𝒙)
        𝝭[i] = 𝒑₀ᵀ𝗠⁻¹*𝒑*𝜙
        ∂𝝭∂x[i] = 𝒑₀ᵀ∂𝗠⁻¹∂x*𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂x*𝜙 + 𝒑₀ᵀ𝗠⁻¹*𝒑*∂𝜙∂x
        ∂𝝭∂y[i] = 𝒑₀ᵀ∂𝗠⁻¹∂y*𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂y*𝜙 + 𝒑₀ᵀ𝗠⁻¹*𝒑*∂𝜙∂y
        ∂𝝭∂z[i] = 𝒑₀ᵀ∂𝗠⁻¹∂z*𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂z*𝜙 + 𝒑₀ᵀ𝗠⁻¹*𝒑*∂𝜙∂z
    end
    return 𝝭, ∂𝝭∂x, ∂𝝭∂y, ∂𝝭∂z
end

## set shape functions
function set𝝭!(aps::Vector{ReproducingKernel{SNode}})
    for ap in aps
        set𝝭!(ap)
    end
end
function set∇𝝭!(aps::Vector{ReproducingKernel{SNode}})
    for ap in aps
        set∇𝝭!(ap)
    end
end

function set𝝭!(ap::ReproducingKernel{SNode})
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        i = ξ.id
        I = ξ.index[i]
        ξ̂ = Node(ξ)
        𝝭 = get𝝭(ap,ξ̂)
        for j in 1:length(𝓒)
            ξ.𝝭[:∂1][I+j] = 𝝭[j]
        end
    end
end

function set∇𝝭!(ap::ReproducingKernel{SNode})
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        i = ξ.id
        I = ξ.index[i]
        ξ̂ = Node(ξ)
        𝝭,∂𝝭∂x,∂𝝭∂y,∂𝝭∂z = get∇𝝭(ap,ξ̂)
        for j in 1:length(𝓒)
            ξ.𝝭[:∂1][I+j] = 𝝭[j]
            ξ.𝝭[:∂x][I+j] = ∂𝝭∂x[j]
            ξ.𝝭[:∂y][I+j] = ∂𝝭∂y[j]
            ξ.𝝭[:∂z][I+j] = ∂𝝭∂z[j]
        end
    end
end

## convert
function ReproducingKernel{𝝃,𝒑,𝑠,𝜙,T}(a::ReproducingKernel{𝝃},b::ReproducingKernel{𝝃}) where {𝝃<:AbstractNode,𝒑,𝒒,𝑠,𝜙}
    𝓒 = a.𝓒
    𝓖 = b.𝓖
    if 𝝃 == SNode
        for ξ in b.𝓖
            for s in keys(ξ.𝝭)
                ξ.𝝭[s] = a.𝓖[1].𝝭[s]
            end
            i = findfirst(η->(η.ξ,η.η,η.γ) == (ξ.ξ,ξ.η,ξ.γ),a.𝓖)
            if i ≠ nothing
                η = a.𝓖[i]
                ξ.index[ξ.id] = η.index[η.id]
                for s in keys(ξ.𝝭)
                    ξ.𝝭[s] = η.𝝭[s]
                end
            else
                n = length(𝓒)
                ξ.index[ξ.id] = haskey(ξ.𝝭,:∂1) ? lastindex(ξ.𝝭[:∂1]) : lastindex(ξ.𝝭[:∂x])
                for s in keys(ξ.𝝭)
                    push!(ξ.𝝭[s],0.0 for i in 1:n)
                end
            end
        end
    end
    𝗠 = a.𝗠
    𝝭 = a.𝝭
    return ReproducingKernel{𝝃,𝒑,𝑠,𝜙,T}(𝓒,𝓖,𝗠,𝝭)
end
