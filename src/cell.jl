
##
function get_shape_functions(ap::T,ξ::Node,g::Val) where T<:AbstractCell
    ξ₁ = ξ[:ξ]
    ξ₂ = ξ[:η]
    ξ₃ = ξ[:γ]
    return get_shape_functions(ap,ξ₁,ξ₂,ξ₃,g)
end
function get_shape_functions(ap::T,ξ::Node,gs::Val...) where T<:AbstractCell
    return (get_shape_functions(ap,ξ,g) for g in gs)
end
function get_coordinates(ap::T,ξ::Node) where T<:AbstractCell
    N = get_shape_functions(ap,ξ,Val(:∂1))
    x = 0.0
    y = 0.0
    z = 0.0
    for i in 1:length(ap.𝓒)
        x += ap.𝓒[i][:x]
        y += ap.𝓒[i][:y]
        z += ap.𝓒[i][:z]
    end
    return x, y, z
end


## AbstractSeg
@inline get_shape_functions(ap::T,ξ::Float64,::Float64,::Float64,g::Val) where T<:AbstractSeg = get_shape_functions(ap,ξ,g)
@inline get_weight(ap::T,ξ::Node) where T<:AbstractSeg = 0.5*ap.L*ξ[:w]

# ---------------- Seg2 -------------------
struct Seg2 <: AbstractSeg
    𝓒::Vector{Node}
    𝓖::Vector{Node}
    L::Float64
end

# constructions of Seg2
function Seg2(𝓒::Vector{Node},𝓖::Vector{Node})
    x₁ = 𝓒[1][:x]
    y₁ = 𝓒[1][:y]
    x₂ = 𝓒[2][:x]
    y₂ = 𝓒[2][:y]
    L = ((x₂-x₁)^2+(y₂-y₁)^2)^0.5
    return Seg2(𝓒,𝓖,L)
end

# actions for Seg2
@inline get_shape_functions(::Seg2,ξ::Float64,::Val{:∂1}) = ((1.0-ξ)*0.5,(1.0+ξ)*0.5)
@inline function get_shape_functions(ap::Seg2,::Float64,::Val{:∂x})
    return (-1.0/ap.L,1.0/ap.L)
end
@inline get_shape_functions(::Seg2,::Float64,::Val{:∂y}) = (0.,0.)
@inline get_shape_functions(::Seg2,::Float64,::Val{:∂z}) = (0.,0.)
