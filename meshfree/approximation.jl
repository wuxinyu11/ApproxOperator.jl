
## Symmetric matrix with packed storge
struct SymMat
    n::Int
    m::Vector{Float64}
end
SymMat(n::Int) = SymMat(n,zeros(Int(n*(n+1)/2)))

function getindex(A::SymMat,i::Int,j::Int)
# @inline function getindex(A::SymMat,i::Int,j::Int)
    i > j ? A.m[Int(j+i*(i-1)/2)] : A.m[Int(i+j*(j-1)/2)]
end

function setindex!(A::SymMat,val::Float64,i::Int,j::Int)
# @inline function setindex!(A::SymMat,val::Float64,i::Int,j::Int)
    A.m[Int(i+j*(j-1)/2)] = val
end
@inline *(A::SymMat,v::AbstractVector{Float64}) = sum(A[1,i]*v[i] for i in 1:A.n)

@inline function -(A::SymMat)
    A.m .= .-A.m
    return A
end

fill!(A::SymMat,val::Float64) = fill!(A.m,val)
@inline fill!(A::SymMat,val::Float64) = fill!(A.m,val)
function inverse!(A::SymMat)
    # n = size(A,1)
    n = A.n
    for i in 1:n
        A[i,i] = 1.0/A[i,i]
        for j in i+1:n
            A[i,j] = - sum(A[i,k]*A[k,j] for k in i:j-1)/A[j,j]
        end
    end
    return A
end

function UUᵀ!(A::SymMat)
    # n = size(A,1)
    n = A.n
    for i in 1:n
        A[i,i] = sum(A[i,k]*A[i,k] for k in i:n)
        for j in i+1:n
            A[i,j] = sum(A[i,k]*A[j,k] for k in j:n)
            # A[j,i] = A[i,j]
        end
    end
    return A
end

function UᵀAU!(A::SymMat,U::SymMat)
    n = A.n
    for i in n:-1:1
        for j in n:-1:i
            A[i,j] = sum(U[k,i]*A[k,l]*U[l,j] for k in 1:i for l in 1:j)
        end
    end
end

function UAUᵀ!(A::SymMat,U::SymMat)
    n = A.n
    for i in 1:n
        for j in i:n
            A[i,j] = sum(U[i,k]*A[k,l]*U[j,l] for k in i:n for l in j:n)
        end
    end
end

function UUᵀAUUᵀ!(A::SymMat,U::SymMat)
    UᵀAU!(A,U)
    UAUᵀ!(A,U)
    return A
end


function cholesky!(A::SymMat)
    n = A.n
    for i in 1:n
        for k in 1:i-1
            A[i,i] -= A[k,i]^2
        end
        A[i,i] = A[i,i]^0.5
        for j in i+1:n
            for k in 1:i-1
                A[i,j] -= A[k,i]A[k,j]
            end
            A[i,j] = A[i,j]/A[i,i]
        end
    end
    return nothing
end

## Spatial Partition
struct RegularGrid<:SpatialPartition
    xmin::Vector{Float64}
    dx::Vector{Float64}
    nx::Vector{Int}
    cells::Vector{Set{Int}}
end

# constructions of RegularGrid
function RegularGrid(x::Vector{PhysicalNode};n::Int=1,γ::Int=1)
    n *= γ
    nₚ  = length(x)
    xmin, xmax = extrema(x[i].x for i in 1:nₚ)
    ymin, ymax = extrema(x[i].y for i in 1:nₚ)
    zmin, zmax = extrema(x[i].z for i in 1:nₚ)
    dx = xmax - xmin
    dy = ymax - ymin
    dz = zmax - zmin
    nd = 0
    pd = 1
    dx > eps() ? (nd += 1;pd *= dx) : dx = 1e-14
    dy > eps() ? (nd += 1;pd *= dy) : dy = 1e-14
    dz > eps() ? (nd += 1;pd *= dz) : dz = 1e-14
    para = (γ*nₚ/pd)^(1/nd)
    nx = ceil(Int, dx * para)
    ny = ceil(Int, dy * para)
    nz = ceil(Int, dz * para)

    cells = Vector{Set{Int}}(undef,nx*ny*nz)
    for i in 1:nx*ny*nz
        cells[i] = Set{Int}()
    end
    for i in 1:nₚ
        ix = floor(Int, (x[i].x - xmin)/dx * nx)
        iy = floor(Int, (x[i].y - ymin)/dy * ny)
        iz = floor(Int, (x[i].z - zmin)/dz * nz)

        ix > nx-1 ? ix = nx-1 : nothing
        iy > ny-1 ? iy = ny-1 : nothing
        iz > nz-1 ? iz = nz-1 : nothing
        for ii in -n:n
            for jj in -n:n
                for kk in -n:n
                    iix = ix + ii
                    iiy = iy + jj
                    iiz = iz + kk

                    iix < 0 ? iix = 0 : nothing
                    iiy < 0 ? iiy = 0 : nothing
                    iiz < 0 ? iiz = 0 : nothing
                    iix > nx-1 ? iix = nx-1 : nothing
                    iiy > ny-1 ? iiy = ny-1 : nothing
                    iiz > nz-1 ? iiz = nz-1 : nothing

                    push!(cells[nx*ny*iiz + nx*iiy + iix + 1], i)
                end
            end
        end
    end
    return RegularGrid([xmin,ymin,zmin],[dx,dy,dz],Int[nx,ny,nz],cells)
end

# actions of RegularGrid
function (rg::RegularGrid)(x::PhysicalNode)
    ix = floor(Int, (x.x[1] - rg.xmin[1])/rg.dx[1] * rg.nx[1])
    iy = floor(Int, (x.x[2] - rg.xmin[2])/rg.dx[2] * rg.nx[2])
    iz = floor(Int, (x.x[3] - rg.xmin[3])/rg.dx[3] * rg.nx[3])

    ix > rg.nx[1]-1 ? ix = rg.nx[1]-1 : nothing
    iy > rg.nx[2]-1 ? iy = rg.nx[2]-1 : nothing
    iz > rg.nx[3]-1 ? iz = rg.nx[3]-1 : nothing
    return rg.cells[rg.nx[1]*rg.nx[2]*iz + rg.nx[1]*iy + ix + 1]
end

function (rg::RegularGrid)(xs::PhysicalNode...)
    indices = Set{Int}()
    for x in xs
        union!(indices,rg(x))
    end
    return indices
end
(rg::RegularGrid)(xs::Vector{Node}) = rg(xs...)

## Basis Function
# ------------ Linear1D ---------------
@inline get_basis_function(::Val{:Linear1D},x::NTuple{3,Float64},::Val{:∂1}) = (1.,x[1])
@inline get_basis_function(::Val{:Linear1D}, ::NTuple{3,Float64},::Val{:∂x}) = (0.,1.)
@inline get_basis_function(::Val{:Linear1D}, ::NTuple{3,Float64},::Val{:∂y}) = (0.,0.)
@inline get_basis_function(::Val{:Linear1D}, ::NTuple{3,Float64},::Val{:∂z}) = (0.,0.)

# ------------ Quadaratic1D ---------------
@inline get_basis_function(::Val{:Quadratic1D},x::NTuple{3,Float64},::Val{:∂1}) = {3,Float64}(1.,x[1],x[1]^2)
@inline get_basis_function(::Val{:Quadratic1D},x::NTuple{3,Float64},::Val{:∂x}) = {3,Float64}(0.,1.,2*x[1])
@inline get_basis_function(::Val{:Quadratic1D}, ::NTuple{3,Float64},::Val{:∂y}) = {3,Float64}(0.,0.,0.)
@inline get_basis_function(::Val{:Quadratic1D}, ::NTuple{3,Float64},::Val{:∂z}) = {3,Float64}(0.,0.,0.)
@inline get_basis_function(::Val{:Quadratic1D}, ::NTuple{3,Float64},::Val{:∂x²}) = {3,Float64}(0.,0.,2.)

# ------------ Cubic1D ---------------
@inline get_basis_function(::Val{:Cubic1D},x::NTuple{3,Float64},::Val{:∂1}) = (1.,x[1],x[1]^2,x[1]^3)
@inline get_basis_function(::Val{:Cubic1D},x::NTuple{3,Float64},::Val{:∂x}) = (0.,1.,2*x[1],3*x[1]^2)
@inline get_basis_function(::Val{:Cubic1D}, ::NTuple{3,Float64},::Val{:∂y}) = (0.,0.,0.,0.)
@inline get_basis_function(::Val{:Cubic1D}, ::NTuple{3,Float64},::Val{:∂z}) = (0.,0.,0.,0.)
@inline get_basis_function(::Val{:Cubic1D},x::NTuple{3,Float64},::Val{:∂x²}) = (0.,0.,2.,6*x[1])

# ------------ Linear2D ---------------
@inline get_basis_function(::Val{:Linear2D},x::NTuple{3,Float64},::Val{:∂1}) = (1.,x[1],x[2])
@inline get_basis_function(::Val{:Linear2D}, ::NTuple{3,Float64},::Val{:∂x}) = (0.,1.,0.)
@inline get_basis_function(::Val{:Linear2D}, ::NTuple{3,Float64},::Val{:∂y}) = (0.,0.,1.)
@inline get_basis_function(::Val{:Linear2D}, ::NTuple{3,Float64},::Val{:∂z}) = (0.,0.,0.)

# ------------ Quadratic2D ---------------
@inline get_basis_function(::Val{:Quadratic2D},x::NTuple{3,Float64},::Val{:∂1}) = (1.,x[1],x[2],x[1]^2,x[1]*x[2],x[2]^2)
@inline get_basis_function(::Val{:Quadratic2D},x::NTuple{3,Float64},::Val{:∂x}) = (0.,1.,0.,2*x[1],x[2],0.)
@inline get_basis_function(::Val{:Quadratic2D},x::NTuple{3,Float64},::Val{:∂y}) = (0.,0.,1.,0.,x[1],2*x[2])
@inline get_basis_function(::Val{:Quadratic2D}, ::NTuple{3,Float64},::Val{:∂z}) = (0.,0.,0.,0.,0.,0.)

# ------------ Cubic2D ---------------
@inline get_basis_function(::Val{:Cubic2D},x::NTuple{3,Float64},::Val{:∂1}) =
(
    1., x[1], x[2], x[1]^2, x[1]*x[2], x[2]^2, x[1]^3, x[1]^2*x[2], x[1]*x[2]^2, x[2]^3
)
@inline get_basis_function(::Val{:Cubic2D},x::NTuple{3,Float64},::Val{:∂x}) =
(
    0., 1., 0., 2*x[1], x[2], 0., 3*x[1]^2, 2*x[1]*x[2], x[2]^2, 0.
)
@inline get_basis_function(::Val{:Cubic2D},x::NTuple{3,Float64},::Val{:∂y}) =
(
    0., 0., 1., 0., x[1], 2*x[2], 0., x[1]^2, 2*x[1]*x[2], 3*x[2]^2
)
@inline get_basis_function(::Val{:Cubic2D},::NTuple{3,Float64},::Val{:∂z}) =
(
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
)

## Kernel Function
# --------------- TensorProductKernel ---------------
struct TensorProductKernel <: KernelFunction
    support_size::NTuple{3,Float64}
    kernel_type::Symbol
end

# constructions of TensorProductKernel
TensorProductKernel(;ss::NTuple{3,Float64}=(1.,1.,1.),kt::Symbol=:CubicSpline) = TensorProductKernel(ss,kt)

# actions of TensorProductKernel
function (kf::TensorProductKernel)(Δx::NTuple{3,Float64},::Val{:∂1})
    sᵢ = kf.support_size
    kt = Val(kf.kernel_type)
    rx = abs(Δx[1])/sᵢ[1]
    ry = abs(Δx[2])/sᵢ[2]
    rz = abs(Δx[3])/sᵢ[3]
    wx = get_kernel(kt,rx,Val(:∂1))
    wy = get_kernel(kt,ry,Val(:∂1))
    wz = get_kernel(kt,rz,Val(:∂1))
    return wx*wy*wz
end

function (kf::TensorProductKernel)(Δx::NTuple{3,Float64},::Val{:∂1},::Val{:∂x})
    sᵢ = kf.support_size
    kt = Val(kf.kernel_type)
    rx = abs(Δx[1])/sᵢ[1]
    ∂rx = sign(Δx[1])/sᵢ[1]
    wx = get_kernel(kt,rx,Val(:∂1))
    ∂wx = get_kernel(kt,rx,Val(:∂r))*∂rx
    return wx, ∂wx
end

function (kf::TensorProductKernel)(Δx::NTuple{3,Float64},::Val{:∂1},::Val{:∂x},::Val{:∂y})
    sᵢ = kf.support_size
    kt = Val(kf.kernel_type)
    rx = abs(Δx[1])/sᵢ[1]
    ry = abs(Δx[2])/sᵢ[2]
    ∂rx = sign(Δx[1])/sᵢ[1]
    ∂ry = sign(Δx[2])/sᵢ[2]
    wx = get_kernel(kt,rx,Val(:∂1))
    wy = get_kernel(kt,ry,Val(:∂1))
    ∂wx = get_kernel(kt,rx,Val(:∂r))*∂rx
    ∂wy = get_kernel(kt,ry,Val(:∂r))*∂ry
    return wx*wy, ∂wx*wy, wx*∂wy
end

function (kf::TensorProductKernel)(Δx::NTuple{3,Float64},::Val{:∂1},::Val{:∂x},::Val{:∂y},::Val{:∂z})
    sᵢ = kf.support_size
    kt = Val(kf.kernel_type)
    rx = abs(Δx[1])/sᵢ[1]
    ry = abs(Δx[2])/sᵢ[2]
    rz = abs(Δx[3])/sᵢ[3]
    ∂rx = sign(Δx[1])/sᵢ[1]
    ∂ry = sign(Δx[2])/sᵢ[2]
    ∂rz = sign(Δx[3])/sᵢ[3]
    wx = get_kernel(kt,rx,Val(:∂1))
    wy = get_kernel(kt,ry,Val(:∂1))
    wz = get_kernel(kt,rz,Val(:∂1))
    ∂wx = get_kernel(kt,rx,Val(:∂r))*∂rx
    ∂wy = get_kernel(kt,ry,Val(:∂r))*∂ry
    ∂wz = get_kernel(kt,rz,Val(:∂r))*∂rz
    return wx*wy*wz, ∂wx*wy*wz, wx*∂wy*wz, wx*wy*∂wz
end

# ----------------- CircularKernel ---------------
struct CircularKernel <: KernelFunction
    support_size::Float64
    kernel_type::Symbol
end

# --------------- Kernel ---------------
function get_kernel(::Val{:CubicSpline},r::Float64,::Val{:∂1})
    if r > 1.0
        return 0.0
    elseif r <= 0.5
        return 2/3 - 4*r^2 +  4*r^3
    else
        return 4/3 - 4*r + 4*r^2 - 4*r^3/3
    end
end

function get_kernel(::Val{:CubicSpline},r::Float64,::Val{:∂r})
    if r > 1.0
        return 0.0
    elseif r <= 0.5
        return - 8*r + 12*r^2
    else
        return - 4   + 8*r - 4*r^2
    end
end

function get_kernel(::Val{:CubicSpline},r::Float64,::Val{:∂r²})
    if r > 1.0
        return 0.0
    elseif r <= 0.5
        return - 8 + 24*r
    else
        return   8 - 8*r
    end
end

## calulate shape functions
function cal_shape_functions(ap::ReproducingKernel,ξ::Union{Float64,AbstractVector{Float64}},::Val{:∂1})
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

function cal_shape_functions(ap::ReproducingKernel,ξ::Union{Float64,AbstractVector{Float64}},::Val{:∂1},::Val{:∂x})
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

function cal_shape_functions(ap::ReproducingKernel,ξ::AbstractVector{Float64},::Val{:∂1},::Val{:∂x},::Val{:∂y})
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

function cal_shape_functions(ap::ReproducingKernel,ξ::AbstractVector{Float64},::Val{:∂1},::Val{:∂x},::Val{:∂y},::Val{:∂z})
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

## MFSpace
mutable struct MFSpace{S<:SpatialPartition,K<:KernelFunction}
    spatialpartition::S
    kernelfunction::K
    𝗠::Dict{Symbol,SymMat}
    𝝭::Dict{Symbol,Vector{Float64}}
end

function (mf::MFSpace)(aps::Vector{T},bf::Val) where T<:Approximator
    for ap in aps
        mf(ap,bf)
    end
end

function (mf::MFSpace)(ap::T,bf::Val) where T<:Approximator
    𝓖 = ap.𝓒
    for ξ in 𝓖
end
