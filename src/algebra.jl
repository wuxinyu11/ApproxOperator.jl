
@inline *(x::NTuple{N,Float64},y::NTuple{N,Float64}) where N = sum(x[i]*y[i] for i in 1:N)
## Symmetric matrix with packed storge
struct SymMat
    n::Int
    m::Vector{Float64}
end
SymMat(n::Int) = SymMat(n,zeros(Int(n*(n+1)/2)))

@inline function getindex(A::SymMat,i::Int,j::Int)
    i > j ? A.m[Int(j+i*(i-1)/2)] : A.m[Int(i+j*(j-1)/2)]
end

@inline function setindex!(A::SymMat,val::Float64,i::Int,j::Int)
    A.m[Int(i+j*(j-1)/2)] = val
end
@inline *(A::SymMat,ReproducingKernel::NTuple{N,Float64}) where N = sum(A[1,i]*v[i] for i in 1:N)
@inline function *(v::NTuple{N,Float64},A::SymMat) where N
    return Tuple(sum(v[i]*A[i,j] for i in 1:N) for j in 1:N)
end

@inline function -(A::SymMat)
    A.m .= .-A.m
    return A
end

@inline fill!(A::SymMat,val::Float64) = fill!(A.m,val)
function inverse!(A::SymMat)
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
# -------------- RegularGrid ------------------
struct RegularGrid<:SpatialPartition
    xmin::Vector{Float64}
    dx::Vector{Float64}
    nx::Vector{Int}
    cells::Vector{Set{Int}}
end

# constructions of RegularGrid
function RegularGrid(x::Vector{Float64},y::Vector{Float64},z::Vector{Float64};n::Int=1,γ::Int=1)
    n *= γ
    nₚ  = length(x)
    xmin, xmax = extrema(x[i] for i in 1:nₚ)
    ymin, ymax = extrema(y[i] for i in 1:nₚ)
    zmin, zmax = extrema(z[i] for i in 1:nₚ)
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
        ix = floor(Int, (x[i] - xmin)/dx * nx)
        iy = floor(Int, (y[i] - ymin)/dy * ny)
        iz = floor(Int, (z[i] - zmin)/dz * nz)

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
function (rg::RegularGrid)(x::Float64,y::Float64,z::Float64)
    ix = floor(Int, (x - rg.xmin[1])/rg.dx[1] * rg.nx[1])
    iy = floor(Int, (y - rg.xmin[2])/rg.dx[2] * rg.nx[2])
    iz = floor(Int, (z - rg.xmin[3])/rg.dx[3] * rg.nx[3])

    ix > rg.nx[1]-1 ? ix = rg.nx[1]-1 : nothing
    iy > rg.nx[2]-1 ? iy = rg.nx[2]-1 : nothing
    iz > rg.nx[3]-1 ? iz = rg.nx[3]-1 : nothing
    return rg.cells[rg.nx[1]*rg.nx[2]*iz + rg.nx[1]*iy + ix + 1]
end

for t in subtypes(SpatialPartition)
    (sp::t)(x::T) where T<:AbstractNode = sp((x.x,x.y,x.z))
    function (sp::t)(xs::T...) where T<:AbstractNode
        indices = Set{Int}()
        for x in xs
            union!(indices,sp(x))
        end
        return indices
    end
    (sp::t)(xs::T) where T<:AbstractVector = sp(xs...)
    function (sp::t)(ap::T) where T<:AbstractElement
        𝓒 = ap.𝓒
        indices = Set{Int}()
        for 𝒙 in 𝓒
            union!(indices,sp(𝒙.x,𝒙.y,𝒙.z))
        end
        union!(𝓒,(Node(i,𝓒[1].data) for i in indices))
    end
    function (sp::t)(aps::Vector{T}) where T<:AbstractElement
        for ap in aps
            sp(ap)
        end
    end
end

## Basis Function
# ------------ Linear1D ---------------
@inline get𝒑(::ReproducingKernel{𝝃,:Linear1D},x::NTuple{3,Float64}) where 𝝃 = (1.,x[1])
@inline get∂𝒑∂x(::ReproducingKernel{𝝃,:Linear1D},::NTuple{3,Float64}) where 𝝃 = (0.,1.)
@inline get∂𝒑∂y(::ReproducingKernel{𝝃,:Linear1D},::Any) where 𝝃 = (0.,0.)
@inline get∂𝒑∂z(::ReproducingKernel{𝝃,:Linear1D},::Any) where 𝝃 = (0.,0.)
@inline get𝒑(::ReproducingKernel{𝝃,:Linear1D},ξ::𝝃) where 𝝃<:AbstractNode = (1.0,0.5*(1.0-ξ.ξ))
@inline get∂𝒑∂x(::ReproducingKernel{𝝃,:Linear1D},::𝝃) where 𝝃<:AbstractNode = (0.0,1.0)

# ------------ Quadaratic1D ---------------
@inline get𝒑(::ReproducingKernel{𝝃,:Quadratic1D},x::NTuple{3,Float64}) where 𝝃 = (1.,x[1],x[1]^2)
@inline get∂𝒑∂x(::ReproducingKernel{𝝃,:Quadratic1D},x::NTuple{3,Float64}) where 𝝃 = (0.,1.,2*x[1])
@inline get∂𝒑∂y(::ReproducingKernel{𝝃,:Quadratic1D},::Any) where 𝝃 = (0.,0.,0.)
@inline get∂𝒑∂z(::ReproducingKernel{𝝃,:Quadratic1D},::Any) where 𝝃 = (0.,0.,0.)
@inline get∂²𝒑∂x²(::ReproducingKernel{𝝃,:Quadratic1D},::Any) where 𝝃 =(0.,0.,2.)
@inline get𝒑(::ReproducingKernel{𝝃,:Quadratic1D},ξ::𝝃) where 𝝃<:AbstractNode = (1.,ξ.ξ,ξ.ξ^2)
@inline get∂𝒑∂x(::ReproducingKernel{𝝃,:Quadratic1D},ξ::𝝃) where 𝝃<:AbstractNode = (0.,-0.5,-ξ.ξ)

# ------------ Cubic1D ---------------
@inline get𝒑(::ReproducingKernel{𝝃,:Cubic1D},x::NTuple{3,Float64}) where 𝝃 = (1.,x[1],x[1]^2,x[1]^3)
@inline get∂𝒑∂x(::ReproducingKernel{𝝃,:Cubic1D},x::NTuple{3,Float64}) where 𝝃 = (0.,1.,2*x[1],3*x[1]^2)
@inline get∂𝒑∂y(::ReproducingKernel{𝝃,:Cubic1D}, ::Any) where 𝝃 = (0.,0.,0.,0.)
@inline get∂𝒑∂z(::ReproducingKernel{𝝃,:Cubic1D}, ::Any) where 𝝃 = (0.,0.,0.,0.)
@inline get∂²𝒑∂x²(::ReproducingKernel{𝝃,:Cubic1D},x::NTuple{3,Float64}) where 𝝃 = (0.,0.,2.,6*x[1])

# ------------ Linear2D ---------------
@inline get𝒑(::ReproducingKernel{𝝃,:Linear2D},x::NTuple{3,Float64}) where 𝝃 = (1.,x[1],x[2])
@inline get∂𝒑∂x(::ReproducingKernel{𝝃,:Linear2D}, ::Any) where 𝝃 = (0.,1.,0.)
@inline get∂𝒑∂y(::ReproducingKernel{𝝃,:Linear2D}, ::Any) where 𝝃 = (0.,0.,1.)
@inline get∂𝒑∂z(::ReproducingKernel{𝝃,:Linear2D}, ::Any) where 𝝃 = (0.,0.,0.)
@inline get𝒑(::ReproducingKernel{𝝃,:Linear2D},ξ::𝝃) where 𝝃<:AbstractNode = (1.,ξ.ξ,ξ.η)

# ------------ Quadratic2D ---------------
@inline get𝒑(::ReproducingKernel{𝝃,:Quadratic2D},x::NTuple{3,Float64}) where 𝝃 = (1.,x[1],x[2],x[1]^2,x[1]*x[2],x[2]^2)
@inline get∂𝒑∂x(::ReproducingKernel{𝝃,:Quadratic2D},x::NTuple{3,Float64}) where 𝝃 = (0.,1.,0.,2*x[1],x[2],0.)
@inline get∂𝒑∂y(::ReproducingKernel{𝝃,:Quadratic2D},x::NTuple{3,Float64}) where 𝝃 = (0.,0.,1.,0.,x[1],2*x[2])
@inline get∂𝒑∂z(::ReproducingKernel{𝝃,:Quadratic2D}, ::Any) where 𝝃 = (0.,0.,0.,0.,0.,0.)
@inline get𝒑(::ReproducingKernel{𝝃,:Quadratic2D},ξ::𝝃) where 𝝃<:AbstractNode = (1.,ξ.ξ,ξ.η,ξ.ξ^2,ξ.ξ*ξ.η,ξ.η^2)
@inline get∂𝒑∂x(::ReproducingKernel{𝝃,:Quadratic2D},ξ::𝝃) where 𝝃<:AbstractNode = (0.,1.,0.,2*ξ.ξ,ξ.η,0.)
@inline get∂𝒑∂y(::ReproducingKernel{𝝃,:Quadratic2D},ξ::𝝃) where 𝝃<:AbstractNode = (0.,0.,1.,0.,ξ.ξ,2*ξ.η)

# ------------ Cubic2D ---------------
@inline get𝒑(::ReproducingKernel{𝝃,:Cubic2D},x::NTuple{3,Float64}) where 𝝃 =
(
    1., x[1], x[2], x[1]^2, x[1]*x[2], x[2]^2, x[1]^3, x[1]^2*x[2], x[1]*x[2]^2, x[2]^3
)
@inline get∂𝒑∂x(::ReproducingKernel{𝝃,:Cubic2D},x::NTuple{3,Float64}) where 𝝃 =
(
    0., 1., 0., 2*x[1], x[2], 0., 3*x[1]^2, 2*x[1]*x[2], x[2]^2, 0.
)
@inline get∂𝒑∂y(::ReproducingKernel{𝝃,:Cubic2D},x::NTuple{3,Float64}) where 𝝃 =
(
    0., 0., 1., 0., x[1], 2*x[2], 0., x[1]^2, 2*x[1]*x[2], 3*x[2]^2
)
@inline get∂𝒑∂z(::ReproducingKernel{𝝃,:Cubic2D},::Any) where 𝝃 =
(
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
)

## Kernel Function
function get𝜙(ap::ReproducingKernel{𝝃,𝒑,:□,𝜙},x::Node,Δx::NTuple{3,Float64}) where {𝝃,𝒑,𝜙}
    rx = abs(Δx[1])/x.s₁
    ry = abs(Δx[2])/x.s₁
    rz = abs(Δx[3])/x.s₁
    wx = get𝜙ᵣ(ap,rx)
    wy = get𝜙ᵣ(ap,ry)
    wz = get𝜙ᵣ(ap,rz)
    return wx*wy*wz
end

function get∂𝜙∂x(ap::ReproducingKernel{𝝃,𝒑,:□,𝜙},x::Node,Δx::NTuple{3,Float64}) where {𝝃,𝒑,𝜙}
    rx = abs(Δx[1])/x.s₁
    ∂rx = sign(Δx[1])/x.s₁
    wx = get𝜙ᵣ(ap,rx)
    ∂wx = get∂𝜙∂r(ap,rx)*∂rx
    return wx, ∂wx
end

function get∇𝜙(::ReproducingKernel{𝝃,𝒑,:□,𝜙},x::Node,Δx::NTuple{3,Float64}) where {𝝃,𝒑,𝜙}
    rx = abs(Δx[1])/x.s₁
    ry = abs(Δx[2])/x.s₂
    rz = abs(Δx[3])/x.s₃
    ∂rx = sign(Δx[1])/x.s₁
    ∂ry = sign(Δx[2])/x.s₂
    ∂rz = sign(Δx[3])/x.s₃
    wx = get𝜙ᵣ(ap,rx)
    wy = get𝜙ᵣ(ap,ry)
    wz = get𝜙ᵣ(ap,rz)
    ∂wx = get∂𝜙∂r(ap,rx)*∂rx
    ∂wy = get∂𝜙∂r(ap,ry)*∂ry
    ∂wz = get∂𝜙∂r(ap,rz)*∂rz
    return wx*wy*wz, ∂wx*wy*wz, wx*∂wy*wz, wx*wy*∂wz
end

## --------------- Kernel ---------------
function get𝜙ᵣ(::ReproducingKernel{𝝃,𝒑,𝑠,:CubicSpline},r::Float64) where {𝝃,𝒑,𝑠}
    if r > 1.0
        return 0.0
    elseif r <= 0.5
        return 2/3 - 4*r^2 +  4*r^3
    else
        return 4/3 - 4*r + 4*r^2 - 4*r^3/3
    end
end

function get∂𝜙∂r(::ReproducingKernel{𝝃,𝒑,𝑠,:CubicSpline},r::Float64) where {𝝃,𝒑,𝑠}
    if r > 1.0
        return 0.0
    elseif r <= 0.5
        return - 8*r + 12*r^2
    else
        return - 4   + 8*r - 4*r^2
    end
end

function get∂²𝜙∂r²(::ReproducingKernel{𝝃,𝒑,𝑠,:CubicSpline},r::Float64) where {𝝃,𝒑,𝑠}
    if r > 1.0
        return 0.0
    elseif r <= 0.5
        return - 8 + 24*r
    else
        return   8 - 8*r
    end
end

## calulate shape functions
function cal𝗠!(ap::ReproducingKernel,x::NTuple{3,Float64})
    𝓒 = ap.𝓒
    𝗠 = ap.𝗠[:∂1]
    n = length(get𝒑(ap.type[1],(0.0,0.0,0.0)))
    fill!(𝗠,0.)
    for xᵢ in 𝓒
        Δx = x - xᵢ
        𝒑 = get𝒑(ap,Δx)
        𝜙 = get𝜙(ap,xᵢ,Δx)
        for I in 1:n
            for J in I:n
                𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
            end
        end
    end
    cholesky!(𝗠)
    U⁻¹ = inverse!(𝗠)
    𝗠⁻¹ = UUᵀ!(U⁻¹)
    return 𝗠⁻¹
end

function cal∂𝗠∂x!(ap::ReproducingKernel,x::NTuple{3,Float64})
    𝓒 = ap.𝓒
    𝗠 = ap.𝗠[:∂1]
    ∂𝗠∂x = ap.𝗠[:∂x]
    n = length(get𝒑(ap.type[1],(0.0,0.0,0.0)))
    fill!(𝗠,0.)
    fill!(∂𝗠∂x,0.)
    for xᵢ in 𝓒
        Δx = x - xᵢ
        𝒑, ∂𝒑∂x = get∂𝒑∂x(ap,Δx)
        𝜙, ∂𝜙∂x = get∂𝜙∂x(ap,xᵢ,Δx)
        for I in 1:n
            for J in I:n
                𝗠[I,J] += w*𝒑[I]*𝒑[J]
                ∂𝗠∂x[I,J] += ∂𝜙∂x*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂x[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂x[J]
            end
        end
    end
    cholesky!(𝗠)
    U⁻¹ = inverse!(𝗠)
    ∂𝗠⁻¹∂x = - UUᵀAUUᵀ!(∂𝗠∂x,U⁻¹)
    𝗠⁻¹ = UUᵀ!(U⁻¹)
    return 𝗠⁻¹, ∂𝗠⁻¹∂x
end

function cal∇𝗠!(ap::ReproducingKernel,x::NTuple{3,Float64})
    𝓒 = ap.𝓒
    𝗠 = ap.𝗠[:∂1]
    ∂𝗠∂x = ap.𝗠[:∂x]
    ∂𝗠∂y = ap.𝗠[:∂y]
    ∂𝗠∂z = ap.𝗠[:∂z]
    n = length(get𝒑(ap.type[1],(0.0,0.0,0.0)))
    fill!(𝗠,0.)
    fill!(∂𝗠∂x,0.)
    fill!(∂𝗠∂y,0.)
    fill!(∂𝗠∂z,0.)
    for xᵢ in 𝓒
        Δx = x - xᵢ
        𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂𝒑∂z = get∇𝒑(ap,Δx)
        𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂𝜙∂z = get∇𝜙(ap,xᵢ,Δx)
        for I in 1:n
            for J in I:n
                𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
                ∂𝗠∂x[I,J] += ∂𝜙∂x*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂x[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂x[J]
                ∂𝗠∂y[I,J] += ∂𝜙∂y*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂y[J]
                ∂𝗠∂z[I,J] += ∂𝜙∂z*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂z[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂z[J]
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

function cal𝗚!(ap::ReproducingKernel)
    𝓖 = ap.𝓖
    𝗚 = ap.𝗠[:∂x]
    n = length(get𝒑(ap,(0.0,0.0,0.0)))
    fill!(𝗚,0.0)
    for ξ in 𝓖
        w = get𝑤(ap,ξ)
        𝒑 = get𝒑(ap,ξ)
        for I in 1:n
            for J in I:n
                𝗚[I,J] += w*𝒑[I]*𝒑[J]
            end
        end
    end
    cholesky!(𝗚)
    U⁻¹ = inverse!(𝗚)
    𝗚⁻¹ = UUᵀ!(U⁻¹)
    return 𝗚⁻¹
end
