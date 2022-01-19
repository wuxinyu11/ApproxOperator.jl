## Counstruction
struct Operator{T,D}
    type::Val{T}
    data::Dict{Symbol,D}
end
Operator(type::Symbol,data::Dict{Symbol,D}) where D = Operator(Val(type),data)
Operator(type::Symbol) = Operator(Val(type))
Operator(type::Symbol,data::DataType) = Operator(Val(type),Dict{Symbol,data}())
Operator(type::Symbol,pair::Pair{Symbol,D}) where D = Operator(Val(type),Dict(pair))

## General Functions
push!(op::Operator,d::Pair{Symbol,D}...) where D<:Any = push!(op.data,d...)

@inline getproperty(op::Operator,f::Symbol) = hasfield(Operator,f) ? getfield(op,f) : getfield(op,:data)[f]
@inline function setproperty!(op::Operator,f::Symbol,x)
    getfield(op,:data)[f] = x
end

@inline function (op::Operator)(aps::Vector{T},k::AbstractMatrix{Float64},f::Vector{Float64}) where T<:Approximator
    for ap in aps
        op(ap,k,f)
    end
end
@inline function (op::Operator)(aps::Vector{T},f::Vector{Float64}) where T<:Approximator
    for ap in aps
        op(ap,f)
    end
end

@inline function (op::Operator)(aps::Vector{T},s::Symbol) where T<:Approximator
    for ap in aps
        op(ap,s)
    end
end
@inline function (op::Operator)(aps::Vector{T}) where T<:Approximator
    for ap in aps
        op(ap)
    end
end

function prescribe!(ap::T,s::Symbol,f::Function) where T<:Approximator
    𝓖 = ap.𝓖
    data = 𝓖[1].data
    if ~haskey(data,s)
        push!(data,s=>similar(data[:w]))
    end
    for ξ in 𝓖
        x = getx(ap,ξ)
        v = f(x...)
        setproperty!(ξ,s,v)
    end
end

function prescribe!(aps::Vector{T},s::Symbol,f::Function) where T<:Approximator
    for ap in aps
        prescribe!(ap,s,f)
    end
end

## Potential Problem
function (op::Operator{:∫∇v∇uvbdΩ})(ap::Approximator,k::AbstractMatrix{Float64},f::AbstractVector{Float64})
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N,B₁,B₂,B₃ = get∇𝝭(ap,ξ)
        w = getw(ap,ξ)
        for i in 1:length(𝓒)
            I = 𝓒[i].id
            for j in 1:length(𝓒)
                J = 𝓒[j].id
                k[I,J] += op.k*(B₁[i]*B₁[j] + B₂[i]*B₂[j] + B₃[i]*B₃[j])*w
            end
            f[I] += N[i]*ξ.b*w
        end
    end
end

function (op::Operator{:∫∇v∇udΩ})(ap::Approximator,k::AbstractMatrix{Float64},f::AbstractVector{Float64})
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        ~,B₁,B₂,B₃ = get∇𝝭(ap,ξ)
        w = getw(ap,ξ)
        for i in 1:length(𝓒)
            I = 𝓒[i].id
            for j in 1:length(𝓒)
                J = 𝓒[j].id
                k[I,J] += op.k*(B₁[i]*B₁[j] + B₂[i]*B₂[j] + B₃[i]*B₃[j])*w
            end
        end
    end
end

function (op::Operator{:∫vbdΩ})(ap::Approximator,k::AbstractMatrix{Float64},f::AbstractVector{Float64})
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        N = get𝝭(ap,ξ)
        w = getw(ap,ξ)
        for i in 1:length(𝓒)
            I = 𝓒[i].id
            f[I] += N[i]*ξ.b*w
        end
    end
end

function (op::Operator{:∫vtdΓ})(ap::Approximator,f::AbstractVector{Float64})
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        w = getw(ap,ξ)
        N = get𝝭(ap,ξ)
        for i in 1:length(𝓒)
            I = 𝓒[i].id
            f[I] += N[i]*ξ.t*w
        end
    end
end

function (op::Operator{:∫vgdΓ})(ap::Approximator,k::AbstractMatrix{Float64},f::AbstractVector{Float64})
    𝓒 = ap.𝓒; 𝓖 = ap.𝓖
    for ξ in 𝓖
        w = getw(ap,ξ)
        N = get𝝭(ap,ξ)
        for i in 1:length(𝓒)
            I = 𝓒[i].id
            for j in 1:length(𝓒)
                J = 𝓒[j].id
                k[I,J] += op.α*N[i]*N[j]*w
            end
            f[I] += op.α*N[i]*ξ.g*w
        end
    end
end

function (op::Operator{:g})(ap::Poi1,k::AbstractMatrix{Float64},f::AbstractVector{Float64};dof::Symbol=:d)
    x = ap.𝓒[1]
    j = x.id
    g = getproperty(x,dof)
    for i in 1:length(f)
        f[i] -= k[i,j]*g
    end
    k[j,:] .= 0.
    k[:,j] .= 0.
    k[j,j] = 1.
    f[j] = g
end

## Meshfree
function (op::Operator{:𝝭})(ap::ReproducingKernel{SNode})
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        i = ξ.id
        push!(ξ.index,ξ.index[i]+length(𝓒))
        ξ̂ = Node(ξ)
        𝝭 = get𝝭(ap,ξ̂)
        push!(ξ.𝝭[:∂1],(𝝭[i] for i in 1:length(𝓒))...)
    end
end
function (op::Operator{:∇𝝭})(ap::ReproducingKernel{SNode})
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        i = ξ.id
        push!(ξ.index,ξ.index[i]+length(𝓒))
        ξ̂ = Node(ξ)
        𝝭,∂𝝭∂x,∂𝝭∂y,∂𝝭∂z = get∇𝝭(ap,ξ̂)
        push!(ξ.𝝭[:∂1],(𝝭[i] for i in 1:length(𝓒))...)
        push!(ξ.𝝭[:∂x],(∂𝝭∂x[i] for i in 1:length(𝓒))...)
        push!(ξ.𝝭[:∂y],(∂𝝭∂y[i] for i in 1:length(𝓒))...)
        push!(ξ.𝝭[:∂z],(∂𝝭∂z[i] for i in 1:length(𝓒))...)
    end
end

function (op::Operator{:𝝭checkrepeat})(ap::ReproducingKernel{SNode})
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        x = getx(ap,ξ)
        i = ξ.id
        if haskey(op.id,x)
            i = op.id[x]
            ids = op.ids[op.index[i]+1:op.index[i+1]]
            index = (findfirst(x->x==ξ_.id,ids) for ξ_ in 𝓒)
            push!(op.index,last(op.index))
        else
            push!(op.id,x=>i)
            push!(op.ids,(ξ_.id for ξ_ in 𝓒)...)
            push!(op.index,last(op.index)+length(𝓒))
            index = 1:length(𝓒)
        end
        push!(ξ.index,last(ξ.index)+length(𝓒))
        ξ̂ = Node(ξ)
        𝝭 = get𝝭(ap,ξ̂)
        push!(ξ.𝝭[:∂1],(i≠nothing ? 𝝭[i] : 0.0 for i in index)...)
    end
end
