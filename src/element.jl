"""
Element{T}<:AbstractElement{T}
"""
struct Element{T}<:AbstractElement{T}
    𝓒::Vector{Node}
    𝓖::Vector{Node}
end

function set𝒙!(aps::Vector{T}) where T<:AbstractElement
    nᵢ = getnᵢ(aps)
    data = getfield(aps[end].𝓖[end],:data)
    push!(data,:x=>(2,zeros(nᵢ)),:y=>(2,zeros(nᵢ)),:z=>(2,zeros(nᵢ)))
    set𝒙!.(aps)
end
function set𝒙!(ap::T) where T<:AbstractElement
    𝓖 = ap.𝓖
    for ξ in 𝓖
        x,y,z = get𝒙(ap,ξ)
        ξ.x = x
        ξ.y = y
        ξ.z = z
    end
end


function set𝑤!(aps::Vector{T}) where T<:AbstractElement
    nᵢ = getnᵢ(aps)
    data = getfield(aps[end].𝓖[end],:data)
    push!(data,:𝑤=>(:𝐺,zeros(nᵢ)))
    set𝑤!.(aps)
end
function set𝑤!(ap::T) where T<:AbstractElement
    𝓖 = ap.𝓖
    for ξ in 𝓖
        𝑤 = get𝑤(ap,ξ)
        ξ.𝑤 = 𝑤
    end
end

"""
get𝐿,get𝐴,get𝑉
"""
function set𝐿!(aps::Vector{T}) where T<:AbstractElement
    nₑ = length(aps)
    push!(getfield(aps[1].𝓖[1],:data),:𝐿=>(3,zeros(nₑ)))
    set𝐿!.(aps)
end
function set𝐿!(ap::T) where T<:AbstractElement
    𝐿 = get𝐿(ap)
    ap.𝓖[1].𝐿 = 𝐿
end
function set𝐴!(aps::Vector{T}) where T<:AbstractElement
    nₑ = length(aps)
    push!(getfield(aps[1].𝓖[1],:data),:𝐴=>(3,zeros(nₑ)))
    set𝐴!.(aps)
end
function set𝐴!(ap::T) where T<:AbstractElement
    𝐴 = get𝐴(ap)
    ap.𝓖[1].𝐴 = 𝐴
end
function set𝑉!(aps::Vector{T}) where T<:AbstractElement
    nₑ = length(aps)
    push!(getfield(aps[1].𝓖[1],:data),:𝑉=>(3,zeros(nₑ)))
    set𝑉!.(aps)
end
function set𝒙ₘ!(aps::Vector{T}) where T<:AbstractElement
    nₑ = length(aps)
    push!(getfield(aps[1].𝓖[1],:data),:xₘ=>(3,zeros(nₑ)),:yₘ=>(3,zeros(nₑ)))
    set𝒙ₘ!.(aps)
end
function set𝒙ₘ!(ap::T) where T<:AbstractElement
    xₘ,yₘ = get𝒙ₘ(ap)
    ap.𝓖[1].xₘ = xₘ
    ap.𝓖[1].yₘ = yₘ
end
function setm2!(aps::Vector{T}) where T<:AbstractElement
    nₑ = length(aps)
    push!(getfield(aps[1].𝓖[1],:data),:m₂₀=>(3,zeros(nₑ)),:m₁₁=>(3,zeros(nₑ)),:m₀₂=>(3,zeros(nₑ)))
    setm2!.(aps)
end
function setm2!(ap::T) where T<:AbstractElement
    m₂₀,m₁₁,m₀₂ = getm2(ap)
    ap.𝓖[1].m₂₀ = m₂₀
    ap.𝓖[1].m₁₁ = m₁₁
    ap.𝓖[1].m₀₂ = m₀₂
end

function get𝐴(ap::T) where T<:AbstractElement{:Vor2}
    𝓒 = ap.𝓒
    nᵥ = length(𝓒)
    𝐴 = 0.0
    for i in 1:nᵥ-1
        x₁ = 𝓒[i].x
        y₁ = 𝓒[i].y
        x₂ = 𝓒[i+1].x
        y₂ = 𝓒[i+1].y
        𝐴 += x₁*y₂-x₂*y₁
    end
    x₁ = 𝓒[nᵥ].x
    y₁ = 𝓒[nᵥ].y
    x₂ = 𝓒[1].x
    y₂ = 𝓒[1].y
    𝐴 += x₁*y₂-x₂*y₁
    𝐴 *= 0.5
    return 𝐴
end

function get𝒙ₘ(ap::AbstractElement)
    𝓒 = ap.𝓒
    nᵥ = length(𝓒)
    𝐴 = ap.𝓖[1].𝐴
    xₘ = 0.0
    yₘ = 0.0
    for i in 1:nᵥ-1
        x₁ = 𝓒[i].x-x₀
        y₁ = 𝓒[i].y-y₀
        x₂ = 𝓒[i+1].x-x₀
        y₂ = 𝓒[i+1].y-y₀
        xₘ += (x₁*y₂-x₂*y₁)*(x₁+x₂)
        yₘ += (x₁*y₂-x₂*y₁)*(y₁+y₂)
    end
    x₁ = 𝓒[nᵥ].x
    y₁ = 𝓒[nᵥ].y
    x₂ = 𝓒[1].x
    y₂ = 𝓒[1].y
    xₘ += (x₁*y₂-x₂*y₁)*(x₁+x₂)
    yₘ += (x₁*y₂-x₂*y₁)*(y₁+y₂)
    xₘ /= 6.0*𝐴
    yₘ /= 6.0*𝐴
    return xₘ,yₘ
end
function getm2(ap::AbstractElement)
    𝓒 = ap.𝓒
    nᵥ = length(𝓒)
    𝐴 = ap.𝓖[1].𝐴
    xₘ = ap.𝓖[1].xₘ
    yₘ = ap.𝓖[1].yₘ
    m₂₀ = 0.0
    m₁₁ = 0.0
    m₀₂ = 0.0
    for i in 1:nᵥ-1
        x₁ = 𝓒[i].x
        y₁ = 𝓒[i].y
        x₂ = 𝓒[i+1].x
        y₂ = 𝓒[i+1].y
        m₂₀ += (x₁*y₂-x₂*y₁)*(x₁^2+x₁*x₂+x₂^2)
        m₁₁ += (x₁*y₂-x₂*y₁)*(2*x₁*y₁+x₁*y₂+x₂*y₁+2*x₂*y₂)
        m₀₂ += (x₁*y₂-x₂*y₁)*(y₁^2+y₁*y₂+y₂^2)
    end
    x₁ = 𝓒[nᵥ].x
    y₁ = 𝓒[nᵥ].y
    x₂ = 𝓒[1].x
    y₂ = 𝓒[1].y
    m₂₀ += (x₁*y₂-x₂*y₁)*(x₁^2+x₁*x₂+x₂^2)
    m₁₁ += (x₁*y₂-x₂*y₁)*(2*x₁*y₁+x₁*y₂+x₂*y₁+2*x₂*y₂)
    m₀₂ += (x₁*y₂-x₂*y₁)*(y₁^2+y₁*y₂+y₂^2)
    m₂₀ /= 12.0*𝐴
    m₁₁ /= 24.0*𝐴
    m₀₂ /= 12.0*𝐴
    return m₂₀-xₘ^2,m₁₁-xₘ*yₘ,m₀₂-yₘ^2
end


"""
setgeometry!(ap::T) where T<:AbstractElement
"""
function setgeometry!(aps::Vector{T}) where T<:AbstractElement
    set𝒙!(aps)
    set𝑤!(aps)
    if T<:AbstractElement{:Seg2}
        set𝐿!(aps)
    elseif T<:AbstractElement{:Tri3}
        set𝐴!(aps)
    elseif T<:AbstractElement{:Tet4}
        set𝑉!(aps)
    end
end

"""
set𝝭!
"""
function set𝝭!(aps::Vector{T}) where T<:AbstractElement
    for ap in aps
        set𝝭!(ap)
    end
end

function set𝝭!(ap::Element{S}) where S
    𝓖 = ap.𝓖
    for ξ in 𝓖
        N = get𝝭(ap,ξ)
        for i in 1:length(ap.𝓒)
            𝝭 = ξ[:𝝭]
            𝝭[i] = N[i]
        end
    end
end

"""
get𝝭(ap::Element,ξ::Node)
"""
# ------------- Poi1 ---------------

for set𝝭 in (:set𝝭!,:set∇𝝭!,:set∇̄𝝭!,:set𝝭̄!,:set∇̃𝝭!)
    @eval begin
        function $set𝝭(aps::Vector{T}) where T<:AbstractElement
            for ap in aps
                𝓖 = ap.𝓖
                for 𝒙 in 𝓖
                    $set𝝭(ap,𝒙)
                end
            end
        end
    end
end

"""
⊆,∩
"""
function issubset(a::T,b::S) where {T<:AbstractElement{:Poi1},S<:AbstractElement{:Seg2}}
    i = findfirst(x->x==a.𝓒[1],b.𝓒)
    return i ≠ nothing && i ≤ 2
end

@inline intersect(a::T,b::T) where T<:AbstractElement = a.𝓒 == b.𝓒 ? a : nothing
@inline function intersect(a::T,b::S) where {T<:AbstractElement{:Seg2},S<:AbstractElement{:Poi1}}
    i = findfirst(x->x==b.𝓒[1],a.𝓒)
    return i ≠ nothing && i ≤ 2 ? a : nothing
end
@inline function intersect(a::T,b::S) where {T<:AbstractElement{:Tri3},S<:AbstractElement{:Poi1}}
    i = findfirst(x->x==b.𝓒[1],a.𝓒)
    return i ≠ nothing && i ≤ 3 ? a : nothing
end
@inline function intersect(a::T,b::S) where {T<:AbstractElement{:Tri3},S<:AbstractElement{:Seg2}}
    i = findfirst(x->x==b.𝓒[1],a.𝓒)
    j = findfirst(x->x==b.𝓒[2],a.𝓒)
    return i ≠ nothing && j ≠ nothing && i ≤ 3 && j ≤ 3 ? a : nothing
end
function intersect(as::Vector{T},bs::Vector{S}) where {T<:AbstractElement,S<:AbstractElement}
    aps = T[]
    for b in bs
        for a in as
            ap = a∩b
            ap ≠ nothing ? push!(aps,ap) : nothing
        end
    end
    return aps
end

"""
getnₚ,getnᵢ,getnₛ
"""
getnₚ(ap::T) where T<:AbstractElement = length(getfield(ap.𝓒[1],:data)[:x][2])
@inline getnₚ(aps::Vector{T}) where T<:AbstractElement = getnₚ(aps[1])

function getnᵢ(aps::Vector{T}) where T<:AbstractElement
    ap = aps[end]
    ξ = ap.𝓖[end]
    return ξ.𝐺
end

function getnₛ(aps::Vector{T}) where T<:AbstractElement
    ap = aps[end]
    ξ = ap.𝓖[end]
    return ξ.𝑠 + length(ap.𝓒)
end

