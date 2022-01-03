## shape function
struct RKShape<:ShapeFunction
    𝝭::Dict{Symbol,Dict{AbstractVector,SparseVector}}
end
## Meshfree
# -------------- PoiM ---------------
struct PoiM{B<:BasisFunction,K<:KernelFunction} <: AbstractPoi
    nodes::Vector{PhysicalNode}
    id::Vector{Int}
    qw::Vector{ParametricNode}
    bf::B
    kf::K
end

# constructions of PoiM
function PoiM(nodes::Vector{Node},id::Vector{Int};bf::BasisFunction=Linear1D(),kf::KernelFunction=TensorProductKernel())
    return PoiM(nodes,id,bf=bf,kf=kf)
end
function PoiM(nodes::Vector{Node},id::Int;bf::BasisFunction=Linear1D(),kf::KernelFunction=TensorProductKernel(),sp::Union{SpatialPartition,Nothing}=nothing)
    if sp ≠ nothing
        id = union!([id],collect(sp(nodes[id])))
    end
    qw = QuadratureRule[:PoiGI1]
    return PoiM(nodes,id,qw,bf,kf)
end

# -------------- SegM ---------------
struct SegM{B<:BasisFunction,K<:KernelFunction} <: AbstractSeg
    nodes :: Vector{PhysicalNode}
    id :: Vector{Int}
    qw::Vector{ParametricNode}
    norm::Float64
    bf::B
    kf::K
end
function SegM(nodes::Vector{PhysicalNode},ids::Vector{Vector{Int}};qw::Symbol=:SegGI2,bf::BasisFunction=Linear1D(),kf::KernelFunction=TensorProductKernel(),sp::Union{SpatialPartition,Nothing}=nothing)
    return [SegM(nodes,id,qw=qw,bf=bf,kf=kf,sp=sp) for id in ids]
end
function SegM(nodes::Vector{PhysicalNode},id::Vector{Int};qw::Symbol=:SegGI2,bf::BasisFunction=Linear1D(),kf::KernelFunction=TensorProductKernel(),sp::Union{SpatialPartition,Nothing}=nothing)
    if sp ≠ nothing
        id = union!(id,collect(sp(nodes[id])))
    end
    L = norm(nodes[id[2]] - nodes[id[1]])
    qw = QuadratureRule[qw]
    return SegM(nodes,id,qw,L,bf,kf)
end

# --------------- TriM ---------------
struct TriM{B<:BasisFunction,K<:KernelFunction} <: AbstractTri
    nodes :: Vector{PhysicalNode}
    id :: Vector{Int}
    qw::Vector{ParametricNode}
    norm :: Float64
    bf:: B
    kf:: K
end

# constructions
function TriM(x::Vector{PhysicalNode},ids::Vector{Vector{Int}};qw::Symbol=:TriGI3,bf::BasisFunction=Linear2D(),kf::KernelFunction=TensorProductKernel(),sp::Union{SpatialPartition,Nothing}=nothing)
    return [TriM(x,id,qw=qw,bf=bf,kf=kf,sp=sp) for id in ids]
end
function TriM(x::Vector{PhysicalNode},id::Vector{Int};qw::Symbol=:TriGI3,bf::BasisFunction=Linear2D(),kf::KernelFunction=TensorProductKernel(),sp::Union{SpatialPartition,Nothing}=nothing)
    if sp ≠ nothing
        id = union!(id,collect(sp(x[id])))
    end
    x1 = x[id[1]].x[1]
    y1 = x[id[1]].x[2]
    z1 = x[id[1]].x[3]
    x2 = x[id[2]].x[1]
    y2 = x[id[2]].x[2]
    z2 = x[id[2]].x[3]
    x3 = x[id[3]].x[1]
    y3 = x[id[3]].x[2]
    z3 = x[id[3]].x[3]
    Ax = 0.5*(y1*z2+y2*z3+y3*z1-y2*z1-y3*z2-y1*z3)
    Ay = 0.5*(z1*x2+z2*x3+z3*x1-z2*x1-z3*x2-z1*x3)
    Az = 0.5*(x1*y2+x2*y3+x3*y1-x2*y1-x3*y2-x1*y3)
    A = (Ax^2 + Ay^2 + Az^2)^0.5
    qw = QuadratureRule[qw]
    return TriM(x,id,qw,A,bf,kf)
end
# -------------- ReproducingKernel ---------------
# actions of ReproducingKernel
ReproducingKernel = Union{SegM{B,K},PoiM{B,K},TriM{B,K}} where {B,K}
function get_shape_functions(ap::ReproducingKernel,ξ::Union{Float64,AbstractVector{Float64}},::Val{:∂1})
    x = get_coordinates(ap,ξ)
    p₀ᵀ𝗠⁻¹ = cal_moment_matrix!(ap,x,Val(:∂1))
    𝝭 = get_shape_function(ap,:∂1)
    for i in 1:get_number_of_indices(ap)
        xᵢ = get_local_node(ap,i)
        Δx = x - xᵢ
        p = get_basis_function(ap.bf,Δx,Val(:∂1))
        w = get_kernel_function(ap.kf,Δx,Val(:∂1))
        𝝭[i] = p₀ᵀ𝗠⁻¹*p*w
    end
    return 𝝭
end

function get_shape_functions(ap::ReproducingKernel,ξ::Union{Float64,AbstractVector{Float64}},::Val{:∂1},::Val{:∂x})
    x = get_coordinates(ap,ξ)
    p₀ᵀ𝗠⁻¹, p₀ᵀ∂𝗠⁻¹∂x = cal_moment_matrix!(ap,x,Val(:∂1),Val(:∂x))
    # 𝝭, ∂𝝭∂x = get_shape_function(ap,:∂1,:∂x)
    𝝭 = get_shape_function(ap,:∂1)
    ∂𝝭∂x = get_shape_function(ap,:∂x)
    for i in 1:get_number_of_indices(ap)
        xᵢ = get_local_node(ap,i)
        Δx = x - xᵢ
        p = get_basis_function(ap.bf,Δx,Val(:∂1))
        ∂p∂x = get_basis_function(ap.bf,Δx,Val(:∂x))
        w, ∂w∂x = get_kernel_function(ap.kf,Δx,Val(:∂1),Val(:∂x))
        𝝭[i] = p₀ᵀ𝗠⁻¹*p*w
        ∂𝝭∂x[i] = p₀ᵀ∂𝗠⁻¹∂x*p*w + p₀ᵀ𝗠⁻¹*∂p∂x*w + p₀ᵀ𝗠⁻¹*p*∂w∂x
    end
    return 𝝭, ∂𝝭∂x
end

function get_shape_functions(ap::ReproducingKernel,ξ::AbstractVector{Float64},::Val{:∂1},::Val{:∂x},::Val{:∂y})
    x = get_coordinates(ap,ξ)
    p₀ᵀ𝗠⁻¹, p₀ᵀ∂𝗠⁻¹∂x, p₀ᵀ∂𝗠⁻¹∂y = cal_moment_matrix!(ap,x,Val(:∂1),Val(:∂x),Val(:∂y))
    # 𝝭, ∂𝝭∂x, ∂𝝭∂y = get_shape_function(ap,:∂1,:∂x,:∂y)
    𝝭 = get_shape_function(ap,:∂1)
    ∂𝝭∂x = get_shape_function(ap,:∂x)
    ∂𝝭∂y = get_shape_function(ap,:∂y)
    for i in 1:get_number_of_indices(ap)
        xᵢ = get_local_node(ap,i)
        Δx = x - xᵢ
        # p, ∂p∂x, ∂p∂y = get_basis_function(ap,Δx,Val(:∂1),Val(:∂x),Val(:∂y))
        p = get_basis_function(ap.bf,Δx,Val(:∂1))
        ∂p∂x = get_basis_function(ap.bf,Δx,Val(:∂x))
        ∂p∂y = get_basis_function(ap.bf,Δx,Val(:∂y))
        w, ∂w∂x, ∂w∂y = get_kernel_function(ap.kf,Δx,Val(:∂1),Val(:∂x),Val(:∂y))
        𝝭[i] = p₀ᵀ𝗠⁻¹*p*w
        ∂𝝭∂x[i] = p₀ᵀ∂𝗠⁻¹∂x*p*w + p₀ᵀ𝗠⁻¹*∂p∂x*w + p₀ᵀ𝗠⁻¹*p*∂w∂x
        ∂𝝭∂y[i] = p₀ᵀ∂𝗠⁻¹∂y*p*w + p₀ᵀ𝗠⁻¹*∂p∂y*w + p₀ᵀ𝗠⁻¹*p*∂w∂y
    end
    return 𝝭, ∂𝝭∂x, ∂𝝭∂y
end

function get_shape_functions(ap::ReproducingKernel,ξ::AbstractVector{Float64},::Val{:∂1},::Val{:∂x},::Val{:∂y},::Val{:∂z})
    x = get_coordinates(ap,ξ)
    p₀ᵀ𝗠⁻¹, p₀ᵀ∂𝗠⁻¹∂x, p₀ᵀ∂𝗠⁻¹∂y, p₀ᵀ∂𝗠⁻¹∂z = cal_moment_matrix!(ap,x,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
    # 𝝭, ∂𝝭∂x, ∂𝝭∂y, ∂𝝭∂z = get_shape_function(ap,:∂1,:∂x,:∂y,:∂z)
    𝝭 = get_shape_function(ap,:∂1)
    ∂𝝭∂x = get_shape_function(ap,:∂x)
    ∂𝝭∂y = get_shape_function(ap,:∂y)
    ∂𝝭∂z = get_shape_function(ap,:∂z)
    for i in 1:get_number_of_indices(ap)
        xᵢ = get_local_node(ap,i)
        Δx = x - xᵢ
        # p, ∂p∂x, ∂p∂y, ∂p∂z = get_basis_function(ap,Δx,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
        p = get_basis_function(ap.bf,Δx,Val(:∂1))
        ∂p∂x = get_basis_function(ap.bf,Δx,Val(:∂x))
        ∂p∂y = get_basis_function(ap.bf,Δx,Val(:∂y))
        ∂p∂z = get_basis_function(ap.bf,Δx,Val(:∂z))
        w, ∂w∂x, ∂w∂y, ∂w∂z = get_kernel_function(ap.kf,Δx,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
        𝝭[i] = p₀ᵀ𝗠⁻¹*p*w
        ∂𝝭∂x[i] = p₀ᵀ∂𝗠⁻¹∂x*p*w + p₀ᵀ𝗠⁻¹*∂p∂x*w + p₀ᵀ𝗠⁻¹*p*∂w∂x
        ∂𝝭∂y[i] = p₀ᵀ∂𝗠⁻¹∂y*p*w + p₀ᵀ𝗠⁻¹*∂p∂y*w + p₀ᵀ𝗠⁻¹*p*∂w∂y
        ∂𝝭∂z[i] = p₀ᵀ∂𝗠⁻¹∂z*p*w + p₀ᵀ𝗠⁻¹*∂p∂z*w + p₀ᵀ𝗠⁻¹*p*∂w∂z
    end
    return 𝝭, ∂𝝭∂x, ∂𝝭∂y, ∂𝝭∂z
end

function cal_moment_matrix!(ap::ReproducingKernel,x::AbstractVector,::Val{:∂1})
    n = get_number_of_basis_function(ap)
    𝗠 = get_moment_matrix(ap,:∂1)
    fill!(𝗠,0.)
    for i in 1:get_number_of_indices(ap)
        xᵢ = get_local_node(ap,i)
        Δx = x - xᵢ
        p = get_basis_function(ap.bf,Δx,Val(:∂1))
        w = get_kernel_function(ap.kf,Δx,Val(:∂1))
        for I in 1:n
            for J in I:n
                𝗠[I,J] += w*p[I]*p[J]
            end
        end
    end
    cholesky!(𝗠)
    U⁻¹ = inverse!(𝗠)
    𝗠⁻¹ = UUᵀ!(U⁻¹)
    return 𝗠⁻¹
end

function cal_moment_matrix!(ap::ReproducingKernel,x::AbstractVector,::Val{:∂1},::Val{:∂x})
    n = get_number_of_basis_function(ap)
    # 𝗠, ∂𝗠∂x = get_moment_matrix(ap,:∂1,:∂x)
    𝗠 = get_moment_matrix(ap,:∂1)
    ∂𝗠∂x = get_moment_matrix(ap,:∂x)
    fill!(𝗠,0.)
    fill!(∂𝗠∂x,0.)
    for i in 1:get_number_of_indices(ap)
        xᵢ = get_local_node(ap,i)
        Δx = x - xᵢ
        # p, ∂p∂x = get_basis_function(ap,Δx,Val(:∂1),Val(:∂x))
        p = get_basis_function(ap.bf,Δx,Val(:∂1))
        ∂p∂x = get_basis_function(ap.bf,Δx,Val(:∂x))
        w, ∂w∂x = get_kernel_function(ap.kf,Δx,Val(:∂1),Val(:∂x))
        for I in 1:n
            for J in I:n
                𝗠[I,J] += w*p[I]*p[J]
                ∂𝗠∂x[I,J] += ∂w∂x*p[I]*p[J] + w*∂p∂x[I]*p[J] + w*p[I]*∂p∂x[J]
            end
        end
    end
    cholesky!(𝗠)
    U⁻¹ = inverse!(𝗠)
    ∂𝗠⁻¹∂x = - UUᵀAUUᵀ!(∂𝗠∂x,U⁻¹)
    𝗠⁻¹ = UUᵀ!(U⁻¹)
    return 𝗠⁻¹, ∂𝗠⁻¹∂x
end

function cal_moment_matrix!(ap::ReproducingKernel,x::AbstractVector,::Val{:∂1},::Val{:∂x},::Val{:∂y},::Val{:∂z})
    n = get_number_of_basis_function(ap)
    # 𝗠, ∂𝗠∂x, ∂𝗠∂y, ∂𝗠∂z = get_moment_matrix(ap,:∂1,:∂x,:∂y,:∂z)
    𝗠 = get_moment_matrix(ap,:∂1)
    ∂𝗠∂x = get_moment_matrix(ap,:∂x)
    ∂𝗠∂y = get_moment_matrix(ap,:∂y)
    ∂𝗠∂z = get_moment_matrix(ap,:∂z)
    fill!(𝗠,0.)
    fill!(∂𝗠∂x,0.)
    fill!(∂𝗠∂y,0.)
    fill!(∂𝗠∂z,0.)
    for i in 1:get_number_of_indices(ap)
        xᵢ = get_local_node(ap,i)
        Δx = x - xᵢ
        # p, ∂p∂x, ∂p∂y, ∂p∂z = get_basis_function(ap,Δx,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
        p = get_basis_function(ap.bf,Δx,Val(:∂1))
        ∂p∂x = get_basis_function(ap.bf,Δx,Val(:∂x))
        ∂p∂y = get_basis_function(ap.bf,Δx,Val(:∂y))
        ∂p∂z = get_basis_function(ap.bf,Δx,Val(:∂z))
        w, ∂w∂x, ∂w∂y, ∂w∂z = get_kernel_function(ap.kf,Δx,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
        for I in 1:n
            for J in I:n
                𝗠[I,J] += w*p[I]*p[J]
                ∂𝗠∂x[I,J] += ∂w∂x*p[I]*p[J] + w*∂p∂x[I]*p[J] + w*p[I]*∂p∂x[J]
                ∂𝗠∂y[I,J] += ∂w∂y*p[I]*p[J] + w*∂p∂y[I]*p[J] + w*p[I]*∂p∂y[J]
                ∂𝗠∂z[I,J] += ∂w∂z*p[I]*p[J] + w*∂p∂z[I]*p[J] + w*p[I]*∂p∂z[J]
            end
        end
    end
    cholesky!(𝗠)
    U⁻¹ = inverse!(𝗠)
    ∂𝗠⁻¹∂x = - UUᵀAUUᵀ!(∂𝗠∂x,U⁻¹)
    ∂𝗠⁻¹∂y = - UUᵀAUUᵀ!(∂𝗠∂y,U⁻¹)
    ∂𝗠⁻¹∂z = - UUᵀAUUᵀ!(∂𝗠∂z,U⁻¹)
    𝗠⁻¹ = UUᵀ!(U⁻¹)
    return 𝗠⁻¹, ∂𝗠⁻¹∂x, ∂𝗠⁻¹∂y, ∂𝗠⁻¹∂z
end

## general functions
# @inline get_basis_function(ap::ReproducingKernel,x::AbstractVector,g::Val) = get_basis_function(ap.bf,x,g)
# @inline get_basis_function(ap::ReproducingKernel,x::AbstractVector,gs::Val...) = (get_basis_function(ap.bf,x,g) for g in gs)
# @inline get_kernel_function(ap::ReproducingKernel,x::AbstractVector,g::Val) = get_kernel_function(ap.kf,x,g)
# @inline get_kernel_function(ap::ReproducingKernel,x::AbstractVector,gs::Val...) = get_kernel_function(ap.kf,x,gs...)
# @inline get_kernel_function(ap::ReproducingKernel,x::AbstractVector,::Val{:∂1}) = get_kernel_function(ap.kf,x,Val(:∂1))
# @inline get_kernel_function(ap::ReproducingKernel,x::AbstractVector,::Val{:∂1},::Val{:∂x}) = get_kernel_function(ap.kf,x,Val(:∂1),Val(:∂x))
# @inline get_kernel_function(ap::ReproducingKernel,x::AbstractVector,::Val{:∂1},::Val{:∂x},::Val{:∂y}) = get_kernel_function(ap.kf,x,Val(:∂1),Val(:∂x),Val(:∂y))
# @inline get_kernel_function(ap::ReproducingKernel,x::AbstractVector,::Val{:∂1},::Val{:∂x},::Val{:∂y},::Val{:∂z}) = get_kernel_function(ap.kf,x,Val(:∂1),Val(:∂x),Val(:∂y),Val(:∂z))
# @inline get_kernel(s::Val,r::Float64,gs::Val...) = (get_kernel(s,r,g) for g in gs)
@inline get_moment_matrix(ap::ReproducingKernel,g::Symbol) = ap.bf.𝗠[g]
# @inline get_moment_matrix(ap::ReproducingKernel,gs::Symbol...) = (ap.bf.𝗠[g] for g in gs)
@inline get_shape_function(ap::ReproducingKernel,g::Symbol) = ap.kf.𝝭[g]
# @inline get_shape_function(ap::ReproducingKernel,gs::Symbol...) = (ap.kf.𝝭[g] for g in gs)
@inline get_number_of_basis_function(ap::ReproducingKernel) = ap.bf.𝗠[:∂1].n
@inline get_number_of_shape_functions(ap::ReproducingKernel) = length(ap.kf.𝝭[:∂1])
