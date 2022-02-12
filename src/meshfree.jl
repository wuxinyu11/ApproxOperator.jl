
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
@inline function setindex!(A::SymMat,val::Float64,i::Int)
    A.m[i] = val
end
@inline *(A::SymMat,v::NTuple{N,Float64}) where N = sum(A[1,i]*v[i] for i in 1:N)
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

## ReproducingKernel
struct ReproducingKernel{𝝃,𝑝,𝑠,𝜙,T}<:AbstractElement{T}
    𝓒::Vector{Node}
    𝓖::Vector{𝝃}
    𝗠::Dict{Symbol,SymMat}
    𝝭::Dict{Symbol,Vector{Float64}}
end

## Basis Function
@inline get∇𝒑(ap::ReproducingKernel,x::Any) = get𝒑(ap,x), get∂𝒑∂x(ap,x), get∂𝒑∂y(ap,x), get∂𝒑∂z(ap,x)
@inline get∇𝒒(ap::ReproducingKernel{𝝃,𝒑,𝑠,𝜙,:Seg2},ξ::Any) where {𝝃<:AbstractNode,𝒑,𝑠,𝜙} = get𝒒(ap,ξ), get∂𝒒∂ξ(ap,ξ)
@inline get∇𝒒(ap::ReproducingKernel{𝝃,𝒑,𝑠,𝜙,:Tri3},ξ::Any) where {𝝃<:AbstractNode,𝒑,𝑠,𝜙} = get𝒒(ap,ξ), get∂𝒒∂ξ(ap,ξ), get∂𝒒∂η(ap,ξ)
@inline get∇𝒒(ap::ReproducingKernel{𝝃,𝒑,𝑠,𝜙,:Tet4},ξ::Any) where {𝝃<:AbstractNode,𝒑,𝑠,𝜙} = get𝒒(ap,ξ), get∂𝒒∂ξ(ap,ξ), get∂𝒒∂η(ap,ξ), get∂𝒒∂γ(ap,ξ)
# ------------ Linear1D ---------------
@inline get𝒑(::ReproducingKernel{𝝃,:Linear1D},x::NTuple{3,Float64}) where 𝝃 = (1.,x[1])
@inline get∂𝒑∂x(::ReproducingKernel{𝝃,:Linear1D},::NTuple{3,Float64}) where 𝝃 = (0.,1.)
@inline get∂𝒑∂y(::ReproducingKernel{𝝃,:Linear1D},::Any) where 𝝃 = (0.,0.)
@inline get∂𝒑∂z(::ReproducingKernel{𝝃,:Linear1D},::Any) where 𝝃 = (0.,0.)
@inline get𝒒(::ReproducingKernel{𝝃,:Linear1D},::Any) where 𝝃<:AbstractNode = (1.0,)
@inline get∂𝒒∂ξ(::ReproducingKernel{𝝃,:Linear1D},::Any) where 𝝃<:AbstractNode = (0.0,)

# ------------ Quadaratic1D ---------------
@inline get𝒑(::ReproducingKernel{𝝃,:Quadratic1D},x::NTuple{3,Float64}) where 𝝃 = (1.,x[1],x[1]^2)
@inline get∂𝒑∂x(::ReproducingKernel{𝝃,:Quadratic1D},x::NTuple{3,Float64}) where 𝝃 = (0.,1.,2*x[1])
@inline get∂𝒑∂y(::ReproducingKernel{𝝃,:Quadratic1D},::Any) where 𝝃 = (0.,0.,0.)
@inline get∂𝒑∂z(::ReproducingKernel{𝝃,:Quadratic1D},::Any) where 𝝃 = (0.,0.,0.)
@inline get∂²𝒑∂x²(::ReproducingKernel{𝝃,:Quadratic1D},::Any) where 𝝃 =(0.,0.,2.)
@inline get𝒒(ap::ReproducingKernel{𝝃,:Quadratic1D},ξ::𝝃) where 𝝃<:AbstractNode = get𝒒(ap,ξ.ξ)
@inline get𝒒(::ReproducingKernel{𝝃,:Quadratic1D},ξ::Float64) where 𝝃<:AbstractNode = (1.0,0.5*(1.0-ξ))
@inline get∂𝒒∂ξ(::ReproducingKernel{𝝃,:Quadratic1D},::Any) where 𝝃<:AbstractNode = (0.0,1.0)

# ------------ Cubic1D ---------------
@inline get𝒑(::ReproducingKernel{𝝃,:Cubic1D},x::NTuple{3,Float64}) where 𝝃 = (1.,x[1],x[1]^2,x[1]^3)
@inline get∂𝒑∂x(::ReproducingKernel{𝝃,:Cubic1D},x::NTuple{3,Float64}) where 𝝃 = (0.,1.,2*x[1],3*x[1]^2)
@inline get∂𝒑∂y(::ReproducingKernel{𝝃,:Cubic1D}, ::Any) where 𝝃 = (0.,0.,0.,0.)
@inline get∂𝒑∂z(::ReproducingKernel{𝝃,:Cubic1D}, ::Any) where 𝝃 = (0.,0.,0.,0.)
@inline get∂²𝒑∂x²(::ReproducingKernel{𝝃,:Cubic1D},x::NTuple{3,Float64}) where 𝝃 = (0.,0.,2.,6*x[1])
@inline get𝒒(ap::ReproducingKernel{𝝃,:Cubic1D},ξ::𝝃) where 𝝃<:AbstractNode = get𝒒(ap,ξ.ξ)
@inline get∂𝒒∂ξ(ap::ReproducingKernel{𝝃,:Cubic1D},ξ::𝝃) where 𝝃<:AbstractNode = get𝒒(ap,ξ.ξ)
@inline get𝒒(::ReproducingKernel{𝝃,:Cubic1D},ξ::Float64) where 𝝃<:AbstractNode = (1.,ξ,ξ^2)
@inline get∂𝒒∂ξ(::ReproducingKernel{𝝃,:Cubic1D},ξ::Float64) where 𝝃<:AbstractNode = (0.,-0.5,-ξ)

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

function get∇𝜙(ap::ReproducingKernel{𝝃,𝒑,:□,𝜙},x::Node,Δx::NTuple{3,Float64}) where {𝝃,𝒑,𝜙}
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
    n = length(get𝒑(ap,(0.0,0.0,0.0)))
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
    n = length(get𝒑(ap,(0.0,0.0,0.0)))
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
    n = length(get𝒑(ap,(0.0,0.0,0.0)))
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
    𝗚 = ap.𝗠[:∇̃]
    n = length(get𝒒(ap,0.0))
    fill!(𝗚,0.0)
    for ξ in 𝓖
        w = get𝑤(ap,ξ)
        𝒒 = get𝒒(ap,ξ)
        for I in 1:n
            for J in I:n
                𝗚[I,J] += w*𝒒[I]*𝒒[J]
            end
        end
    end
    cholesky!(𝗚)
    U⁻¹ = inverse!(𝗚)
    𝗚⁻¹ = UUᵀ!(U⁻¹)
    return 𝗚⁻¹
end

function cal𝗚!(ap::ReproducingKernel{𝝃,𝒑,𝑠,𝜙,:Seg2}) where {𝝃<:AbstractNode,𝒑,𝑠,𝜙}
    𝗚⁻¹ = ap.𝗠[:∇̃]
    fill!(𝗚⁻¹,0.0)
    𝐿 = get𝐿(ap)
    𝗚⁻¹[1] =  4.0/𝐿
    𝗚⁻¹[2] = -6.0/𝐿
    𝗚⁻¹[3] = 12.0/𝐿
    return 𝗚⁻¹
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

function get∇𝑛𝝭(ap::ReproducingKernel{𝝃,𝒑,𝑠,𝜙,:Seg2},ξ::Any) where {𝝃<:AbstractNode,𝒑,𝑠,𝜙}
    N,B₁ = get∇𝝭(ap,ξ)
    n₁ = get𝒏(ap,ξ)
    B = B₁*n₁
    return N, B
end
## set shape functions
function set𝝭!(aps::Vector{T}) where T<:ReproducingKernel{SNode}
    for ap in aps
        set𝝭!(ap)
    end
end
function set∇𝝭!(aps::Vector{T}) where T<:ReproducingKernel{SNode}
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

## shape functions for SNode
function get𝝭(ap::ReproducingKernel,ξ::SNode)
    𝝭 = ap.𝝭[:∂1]
    i = ξ.id
    index = ξ.index
    for j in 1:length(ap.𝓒)
        𝝭[j] = ξ.𝝭[:∂1][index[i]+j]
    end
    return 𝝭
end

function get∂𝝭∂x(ap::ReproducingKernel,ξ::SNode)
    ∂𝝭∂x = ap.𝝭[:∂x]
    i = ξ.id
    index = ξ.index
    for j in 1:length(ap.𝓒)
        ∂𝝭∂x[j] = ξ.𝝭[:∂x][index[i]+j]
    end
    return ∂𝝭∂x
end

function get∇𝝭(ap::ReproducingKernel,ξ::SNode)
    𝝭 = ap.𝝭[:∂1]
    ∂𝝭∂x = ξ.𝝭[:∂x]
    ∂𝝭∂y = ξ.𝝭[:∂y]
    ∂𝝭∂z = ξ.𝝭[:∂z]
    i = ξ.id
    index = ξ.index
    for j in 1:length(ap.𝓒)
        𝝭[j] = ξ.𝝭[:∂1][index[i]+j]
        ∂𝝭∂x[j] = ξ.𝝭[:∂x][index[i]+j]
        ∂𝝭∂y[j] = ξ.𝝭[:∂y][index[i]+j]
        ∂𝝭∂z[j] = ξ.𝝭[:∂z][index[i]+j]
    end
    return 𝝭, ∂𝝭∂x, ∂𝝭∂y, ∂𝝭∂z
end

## RK gradient smoothing
function set∇̃𝝭!(aps::Vector{T}) where T<:ReproducingKernel
    for ap in aps
        set∇̃𝝭!(ap)
    end
end

function set∇̃𝝭!(gps::Vector{T},aps::Vector{S}) where {T<:ReproducingKernel,S<:ReproducingKernel}
    if length(gps) ≠ length(aps)
        error("Miss match element numbers")
    else
        for i in 1:length(gps)
            set∇̃𝝭!(gps[i],aps[i])
        end
    end
end
function set∇̃𝝭!(as::Vector{T},bs::Vector{S},cs::Vector{R}) where {T<:ReproducingKernel,S<:ReproducingKernel,R<:ReproducingKernel}
    if length(as) ≠ length(bs) || length(bs) ≠ length(cs)
        error("Miss match element numbers")
    else
        for i in 1:length(as)
            set∇̃𝝭!(as[i],bs[i],cs[i])
        end
    end
end
set∇̃𝝭!(ap::T) where T<:ReproducingKernel{SNode} = set∇̃𝝭!(ap,ap)
function set∇̃𝝭!(gp::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Seg2},ap::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Seg2}) where {𝒑,𝑠,𝜙}
    n₁ =  1.0
    n₂ = -1.0
    𝗚⁻¹ = cal𝗚!(gp)
    𝓒 = gp.𝓒
    𝓖 = gp.𝓖
    for ξ̂ in 𝓖
        𝒒̂ = get𝒒(gp,ξ̂)
        𝒒̂ᵀ𝗚⁻¹ = 𝒒̂*𝗚⁻¹
        ∂𝝭∂x = gp.𝝭[:∂x]
        fill!(∂𝝭∂x,0.0)
        for ξ in ap.𝓖
            w = ξ.w/2
            wᵇ = ξ.wᵇ
            nᵇ₁ = 0.0
            nᵇ₁ += ξ.ξ ==  1.0 ? n₁ : 0.0
            nᵇ₁ += ξ.ξ == -1.0 ? n₂ : 0.0
            𝝭 = get𝝭(ap,ξ)
            𝒒, ∂𝒒∂ξ = get∇𝒒(gp,ξ)
            W₁ = 𝒒̂ᵀ𝗚⁻¹*𝒒*nᵇ₁*wᵇ + 𝒒̂ᵀ𝗚⁻¹*∂𝒒∂ξ*n₁*w
            for i in 1:length(𝓒)
                ∂𝝭∂x[i] += 𝝭[i]*W₁
            end
        end
        for i in 1:length(𝓒)
            ξ̂.𝝭[:∂x][ξ̂.index[ξ̂.id]+i] = ∂𝝭∂x[i]
        end
    end
end

function set∇̃𝝭!(gp::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Tri3},ap::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Tri3}) where {𝒑,𝑠,𝜙}
    x₁ = gp.𝓒[1].x;y₁ = gp.𝓒[1].y
    x₂ = gp.𝓒[2].x;y₂ = gp.𝓒[2].y
    x₃ = gp.𝓒[3].x;y₃ = gp.𝓒[3].y
    n₁₁ = y₃-y₂;n₂₁ = y₁-y₃;n₃₁ = y₂-y₁
    n₁₂ = x₃-x₂;n₂₂ = x₃-x₁;n₃₂ = x₁-x₂
    𝗚⁻¹ = cal𝗚!(gp)
    𝓒 = gp.𝓒
    𝓖 = gp.𝓖
    for ξ̂ in 𝓖
        𝒒̂ = get𝒒(gp,ξ̂)
        𝒒̂ᵀ𝗚⁻¹ = 𝒒̂*𝗚⁻¹
        ∂𝝭∂x = gp.𝝭[:∂x]
        ∂𝝭∂y = gp.𝝭[:∂y]
        fill!(∂𝝭∂x,0.0)
        fill!(∂𝝭∂y,0.0)
        for ξ in ap.𝓖
            w = ξ.w
            wᵇ = ξ.wᵇ
            n₁ = ξ.n₁
            n₂ = ξ.n₂
            𝝭 = get𝝭(ap,ξ)
            𝒒, ∂𝒒∂ξ, ∂𝒒∂η = get∇𝒒(gp,ξ)
            𝒒̂ᵀ𝗚⁻¹∂𝒒 =  𝒒̂ᵀ𝗚⁻¹*∂𝒒
            𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ = 𝒒̂ᵀ𝗚⁻¹*∂𝒒∂ξ
            𝒒̂ᵀ𝗚⁻¹∂𝒒∂η = 𝒒̂ᵀ𝗚⁻¹*∂𝒒∂η
            nᵇ₁ = 0.0;nᵇ₂ = 0.0
            nᵇ₁ += ξ.ξ == 0.0 ? n₁₁ : 0.0
            nᵇ₁ += ξ.η == 0.0 ? n₂₁ : 0.0
            nᵇ₂ += ξ.ξ == 0.0 ? n₁₂ : 0.0
            nᵇ₂ += ξ.η == 0.0 ? n₂₂ : 0.0
            b₁ = 𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ*n₁₁ + 𝒒̂ᵀ𝗚⁻¹∂𝒒∂η*n₂₁
            b₂ = 𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ*n₁₂ + 𝒒̂ᵀ𝗚⁻¹∂𝒒∂η*n₂₂
            W₁ = 𝒒̂ᵀ𝗚⁻¹𝒒*nᵇ₁*wᵇ + b₁*w/2
            W₂ = 𝒒̂ᵀ𝗚⁻¹𝒒*nᵇ₂*wᵇ + b₂*w/2
            for i in 1:length(𝓒)
                ∂𝝭∂x[i] += 𝝭[i]*W₁
                ∂𝝭∂y[i] += 𝝭[i]*W₂
            end
        end
        for i in 1:length(𝓒)
            ξ̂.𝝭[:∂x][ξ̂.index[ξ̂.id]+i] = ∂𝝭∂x[i]
            ξ̂.𝝭[:∂y][ξ̂.index[ξ̂.id]+i] = ∂𝝭∂y[i]
        end
    end
end

function set∇̃𝝭!(gp::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Tet4},ap::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Tet4}) where {𝒑,𝑠,𝜙}
    𝗚⁻¹ = cal𝗚!(gp)
    𝓒 = gp.𝓒
    𝓖 = gp.𝓖
    for ξ̂ in 𝓖
        𝒒̂ = get𝒒(gp,ξ̂)
        𝒒̂ᵀ𝗚⁻¹ = 𝒒̂*𝗚⁻¹
        ∂𝝭∂x = gp.𝝭[:∂x]
        ∂𝝭∂y = gp.𝝭[:∂y]
        ∂𝝭∂z = gp.𝝭[:∂z]
        fill!(∂𝝭∂x,0.0)
        fill!(∂𝝭∂y,0.0)
        fill!(∂𝝭∂z,0.0)
        for ξ in ap.𝓖
            w = ξ.w
            wᵇ = ξ.wᵇ
            n₁ = ξ.n₁
            n₂ = ξ.n₂
            n₃ = ξ.n₃
            𝝭 = get𝝭(ap,ξ)
            𝒒, ∂𝒒∂ξ, ∂𝒒∂η, ∂𝒒∂γ = get∇𝒒(gp,ξ)
            b = 𝒒̂ᵀ𝗚⁻¹*∂𝒒∂ξ*n₁ + 𝒒̂ᵀ𝗚⁻¹*∂𝒒∂η*n₂ + 𝒒̂ᵀ𝗚⁻¹*∂𝒑∂γ*n₃
            W₁ = 𝒒̂ᵀ𝗚⁻¹*𝒒*n₁*wᵇ + b*w/3
            W₂ = 𝒒̂ᵀ𝗚⁻¹*𝒒*n₂*wᵇ + b*w/3
            W₃ = 𝒒̂ᵀ𝗚⁻¹*𝒒*n₃*wᵇ + b*w/3
            for i in 1:length(𝓒)
                ∂𝝭∂x[i] += 𝝭[i]*W₁
                ∂𝝭∂y[i] += 𝝭[i]*W₂
                ∂𝝭∂z[i] += 𝝭[i]*W₃
            end
        end
        for i in 1:length(𝓒)
            ξ̂.𝝭[:∂x][ξ̂.index[ξ̂.id]+i] = ∂𝝭∂x[i]
            ξ̂.𝝭[:∂y][ξ̂.index[ξ̂.id]+i] = ∂𝝭∂y[i]
            ξ̂.𝝭[:∂z][ξ̂.index[ξ̂.id]+i] = ∂𝝭∂z[i]
        end
    end
end

function setg̃!(gps::Vector{T},aps::Vector{S}) where{T<:ReproducingKernel,S<:ReproducingKernel}
    if length(gps) ≠ length(aps)
        error("Miss match element numbers")
    else
        for i in 1:length(gps)
            setg̃!(gps[i],aps[i])
        end
    end
end

function setg̃!(gp::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Seg2},ap::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Seg2}) where {𝒑,𝑠,𝜙}
    n₁ =  1.0
    n₂ = -1.0
    𝗚⁻¹ = cal𝗚!(gp)
    𝓒 = gp.𝓒
    𝓖 = gp.𝓖
    for ξ̂ in 𝓖
        𝒒̂ = get𝒒(gp,ξ̂)
        𝒒̂ᵀ𝗚⁻¹ = 𝒒̂*𝗚⁻¹
        𝝭 = gp.𝝭[:∂1]
        g̃ = 0.0
        fill!(𝝭,0.0)
        for ξ in ap.𝓖
            w = ξ.w
            n = 0.0
            n += ξ.ξ ==  1.0 ? n₁ : 0.0
            n += ξ.ξ == -1.0 ? n₂ : 0.0
            𝝭 = get𝝭(ap,ξ)
            g = ξ.g
            𝒒 = get𝒒(gp,ξ)
            W₁ = 𝒒̂ᵀ𝗚⁻¹*𝒒*n*w
            for i in 1:length(𝓒)
                𝝭[i] += 𝝭[i]*W₁
            end
            g̃ += 𝒒̂ᵀ𝗚⁻¹*𝒒*g*n*w
        end
        ξ̂.g = g̃
        for i in 1:length(𝓒)
            ξ̂.𝝭[:∂1][ξ̂.index[ξ̂.id]+i] = 𝝭[i]
        end
    end
end

@inline function set∇̃𝝭!(a::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Seg2},b::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Seg2},c::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Seg2}) where {𝒑,𝑠,𝜙}
    set∇̃𝝭!(b,c)
    setg̃!(a,b)
end
## convert
function ReproducingKernel{𝝃,𝒑,𝑠,𝜙,T}(a::ReproducingKernel{𝜼,𝒒}) where {𝝃<:AbstractNode,𝜼<:AbstractNode,𝒑,𝒒,𝑠,𝜙,T}
    𝓒 = a.𝓒
    𝓖 = 𝝃[]
    𝗠 = a.𝗠
    𝝭 = a.𝝭
    b = ReproducingKernel{𝝃,𝒑,𝑠,𝜙,T}(𝓒,𝓖,𝗠,𝝭)
    if 𝒑 ≠ 𝒒
        n = length(get𝒑(b,(0.0,0.0,0.0)))
        for s in keys(𝗠)
            𝗠[s] = SymMat(n)
        end
    end
    return b
end

function ReproducingKernel{𝝃,𝒑,𝑠,𝜙,T}(as::Vector{S}) where {𝝃<:AbstractNode,𝒑,𝑠,𝜙,T,S<:ReproducingKernel}
    aps = ReproducingKernel{𝝃,𝒑,𝑠,𝜙,T}[]
    for a in as
        push!(aps,ReproducingKernel{𝝃,𝒑,𝑠,𝜙,T}(a))
    end
    return aps
end

function ReproducingKernel{𝝃,𝒑,𝑠,𝜙,T}(a::Element{S},𝗠::Dict{Symbol,SymMat},𝝭::Dict{Symbol,Vector{Float64}}) where {𝝃<:AbstractNode,𝜼<:AbstractNode,𝒑,𝒒,𝑠,𝜙,T,S}
    𝓒 = a.𝓒
    𝓖 = 𝝃[]
    return ReproducingKernel{𝝃,𝒑,𝑠,𝜙,T}(𝓒,𝓖,𝗠,𝝭)
end

function ReproducingKernel{𝝃,𝒑,𝑠,𝜙,T}(as::Vector{Element{S}},sp::Union{Nothing,SpatialPartition}=nothing;renumbering::Bool=false) where {𝝃<:AbstractNode,𝒑,𝑠,𝜙,T,S}
    aps = ReproducingKernel{𝝃,𝒑,𝑠,𝜙,T}[]
    𝗠 = Dict{Symbol,SymMat}()
    𝝭 = Dict{Symbol,Vector{Float64}}()
    if renumbering
        index, data = renumber(aps)
        for a in as
            𝓒 = [Node(index[x.id],data) for x in a.𝓒]
            𝓖 = Node[]
            ap = ReproducingKernel{𝝃,𝒑,𝑠,𝜙,T}(𝓒,𝓖,𝗠,𝝭)
            sp ≠ nothing ? sp(ap) : nothing
            push!(aps,ap)
        end
    else
        for a in as
            ap = ReproducingKernel{𝝃,𝒑,𝑠,𝜙,T}(a,𝗠,𝝭)
            sp ≠ nothing ? sp(ap) : nothing
            push!(aps,ap)
        end
    end
    return aps
end
function ReproducingKernel{𝝃,𝒑,𝑠,𝜙,T}(a::A,b::B;sharing::Bool=false) where {𝝃<:AbstractNode,𝒑,𝑠,𝜙,T,A<:ReproducingKernel{𝝃},B<:ReproducingKernel{𝝃}}
    𝓒 = a.𝓒
    𝓖 = get𝓖(a,b)
    if 𝓖 ≠ nothing
        if 𝝃 == SNode
            sharing ? glue(a.𝓖,b.𝓖) : addindex(𝓖,length(a.𝓒)-length(b.𝓒))
        end
        𝗠 = a.𝗠
        𝝭 = a.𝝭
        return ReproducingKernel{𝝃,𝒑,𝑠,𝜙,T}(𝓒,𝓖,𝗠,𝝭)
    else
        return nothing
    end
end

function ReproducingKernel{𝝃,𝒑,𝑠,𝜙,T}(as::Vector{A},bs::Vector{B};sharing::Bool=false) where {𝝃<:AbstractNode,𝒑,𝑠,𝜙,T,A<:ReproducingKernel,B<:ReproducingKernel}
    aps = ReproducingKernel{𝝃,𝒑,𝑠,𝜙,T}[]
    for a in as
        for b in bs
            ap = ReproducingKernel{𝝃,𝒑,𝑠,𝜙,T}(a,b,sharing=sharing)
            ap ≠ nothing ? push!(aps,ap) : nothing
        end
    end
    return aps
end

function addindex(𝓖::Vector{SNode},n::Int)
    nₜ = length(𝓖)*n
    index = 𝓖[1].index
    𝝭 = 𝓖[1].𝝭
    for s in keys(𝝭)
        append!(𝝭[s],zeros(nₜ))
    end
    for ξ in 𝓖
        for i in 1:length(index)-ξ.id
            index[ξ.id+i] += i*n
        end
    end
end

function glue(𝓖₁::Vector{SNode},𝓖₂::Vector{SNode})
    for ξ in 𝓖₂
        for s in keys(ξ.𝝭)
            ξ.𝝭[s] = 𝓖₁[1].𝝭[s]
        end
        i = findfirst(η->(η.ξ,η.η,η.γ) == (ξ.ξ,ξ.η,ξ.γ),𝓖₁)
        if i ≠ nothing
            η = 𝓖₁[i]
            ξ.index[ξ.id] = η.index[η.id]
        else
            η = 𝓖₁[1]
            n = η.index[η.id+1] - η.index[η.id]
            nₜ = 0
            for s in keys(ξ.𝝭)
                nₜ = length(ξ.𝝭[s])
                append!(ξ.𝝭[s],zeros(n))
            end
            ξ.index[ξ.id] = nₜ
        end
    end
end
