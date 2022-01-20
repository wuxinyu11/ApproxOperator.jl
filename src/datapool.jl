
mutable struct Nodes<:DataPool
    nₛ::Int
    nₚ::Int
    nᵢ::Int
    nₑ::Int
    physicaldatas::Dict{Symbol,Float64}
    parametricdatas::Dict{Symbol,Float64}
end

Nodes() = Nodes(0,0,0,0,Dict{Symbol,Float64}(),Dict{Symbol,Float64}())
Nodes(nₛ::Int,nₚ::Int,datas)

function (dp::Nodes)(qtype::Symbol)
    for qr in QuadratureRule[qtype]
        push!(dp.parametricdatas["w"]),qr[1]
        push!(dp.parametricdatas["ξ"]),qr[2]
        if nₛ > 1
            push!(dp.parametricdatas["η"],qr[3])
        end
        if nₛ == 3
            push!(dp.parametricdatas["γ"],qr[4])
        end
        if length(qr) > nₛ+1
            push!(dp.parametricdatas["wᵇ"],last(qr))
        end
    end
    nᵢ = dp.nᵢ
    dp.nᵢ += length(qr)
    return [Node(i,dp.parametricdatas) for i in nᵢ+1:nᵢ+length(qr)]
end

mutable struct SNodes<:DataPool
    nₛ::Int
    nₚ::Int
    nᵢ::Int
    nₑ::Int
    physicaldatas::Dict{Symbol,Float64}
    parametricdatas::Dict{Symbol,Float64}
    𝗠ᵗ::Dict{Symbol,SymMat}
    𝝭ᵗ::Dict{Symbol,Vector{Float64}}
    index::Vector{Int}
    𝝭::Dict{Symbol,Vector{Float64}}
    type::Tuple{Val{𝒑},Val{𝑠},Val{𝜙}}
end
