@inline getdata(ap::T,i::Int,::Val{:∇𝝭}) where T<:Approximator = (ap[:𝝭,i],ap[:∂𝝭∂x,i],ap[:∂𝝭∂y,i],ap[:∂𝝭∂z,i])

## AbstractSeg
(ap::T)(ξ::Node,s::Symbol) where T<:AbstractSeg = ap(ξ,Val(s))
(ap::T)(ξ::Node,::Val{:w}) where T<:AbstractSeg = 0.5*ap.L*ξ[2]
# ----------------- Seg2 -----------------
struct Seg2<:AbstractSeg
    𝓒::Vector{Node{3,Float64}}
    𝓖::Vector{Node{2,Float64}}
    L::Float64
end

# constructions of Seg2
function Seg2(𝓒::Vector{Node{3,Float64}},𝓖::Vector{Node{2,Float64}})
    x₁ = 𝓒[1][1]
    y₁ = 𝓒[1][2]
    x₂ = 𝓒[2][1]
    y₂ = 𝓒[2][2]
    L = ((x₂-x₁)^2+(y₂-y₁)^2)^0.5
    return Seg2(𝓒,𝓖,L)
end

# actions for Seg2
@inline function getdata(ap::Seg2,i::Int,::Val{:𝝭})
    ξ = ap.𝓖[i][1]
    return (0.5*(1-ξ),0.5*(1+ξ))
end
@inline getdata(ap::Seg2,::Int,::Val{:∂𝝭∂x}) = (-1.0/ap.L,1.0/ap.L)
@inline getdata(  ::Seg2,::Int,::Val{:∂𝝭∂y}) = (0.,0.)
@inline getdata(  ::Seg2,::Int,::Val{:∂𝝭∂z}) = (0.,0.)

##
for t in subtypes(Approximator)
