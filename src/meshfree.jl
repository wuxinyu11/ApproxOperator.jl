
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
    for j in 1:N
        A[1,j] = sum(v[i]*A[i,j] for i in 1:N)
    end
    return A
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
function UᵀAU!(B::SymMat,A::SymMat,U::SymMat)
    n = A.n
    for i in n:-1:1
        for j in n:-1:i
            B[i,j] = sum(U[k,i]*A[k,l]*U[l,j] for k in 1:i for l in 1:j)
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

function UUᵀAUUᵀ!(B::SymMat,A::SymMat,U::SymMat)
    UᵀAU!(B,A,U)
    UAUᵀ!(B,U)
    return B
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
    (sp::t)(x::T) where T<:AbstractNode = sp(x.x,x.y,x.z)
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
@inline get∇²𝒑(ap::ReproducingKernel,x::Any) = get𝒑(ap,x), get∂𝒑∂x(ap,x), get∂𝒑∂y(ap,x), get∂²𝒑∂x²(ap,x), get∂²𝒑∂x∂y(ap,x), get∂²𝒑∂y²(ap,x), get∂𝒑∂z(ap,x), get∂²𝒑∂x∂z(ap,x), get∂²𝒑∂y∂z(ap,x), get∂²𝒑∂z²(ap,x)
# @inline get∇²𝒑(ap::ReproducingKernel,x::Any) = get𝒑(ap,x), get∂𝒑∂x(ap,x), get∂𝒑∂y(ap,x), get∂²𝒑∂x²(ap,x), get∂²𝒑∂x∂y(ap,x), get∂²𝒑∂y²(ap,x)
@inline get∇³𝒑(ap::ReproducingKernel,x::Any) = get𝒑(ap,x), get∂𝒑∂x(ap,x), get∂𝒑∂y(ap,x), get∂²𝒑∂x²(ap,x), get∂²𝒑∂x∂y(ap,x), get∂²𝒑∂y²(ap,x), get∂³𝒑∂x³(ap,x), get∂³𝒑∂x²∂y(ap,x), get∂³𝒑∂x∂y²(ap,x), get∂³𝒑∂y³(ap,x)
@inline get∇𝒑₁(ap::ReproducingKernel{𝝃,𝒑,𝑠,𝜙,:Seg2},ξ::Any) where {𝝃<:AbstractNode,𝒑,𝑠,𝜙} = get𝒑₁(ap,ξ), get∂𝒑₁∂ξ(ap,ξ)
@inline get∇𝒑₁(ap::ReproducingKernel{𝝃,𝒑,𝑠,𝜙,:Tri3},ξ::Any) where {𝝃<:AbstractNode,𝒑,𝑠,𝜙} = get𝒑₁(ap,ξ), get∂𝒑₁∂ξ(ap,ξ), get∂𝒑₁∂η(ap,ξ)
@inline get∇²𝒑₂(ap::ReproducingKernel{𝝃,𝒑,𝑠,𝜙,:Tri3},ξ::Any) where {𝝃<:AbstractNode,𝒑,𝑠,𝜙} = get𝒑₂(ap,ξ), get∂𝒑₂∂ξ(ap,ξ), get∂𝒑₂∂η(ap,ξ), get∂²𝒑₂∂ξ²(ap,ξ), get∂²𝒑₂∂ξ∂η(ap,ξ), get∂²𝒑₂∂η²(ap,ξ)
@inline get∇𝒒(ap::ReproducingKernel{𝝃,𝒑,𝑠,𝜙,:Tet4},ξ::Any) where {𝝃<:AbstractNode,𝒑,𝑠,𝜙} = get𝒑₁(ap,ξ), get∂𝒑₁∂ξ(ap,ξ), get∂𝒑₁∂η(ap,ξ), get∂𝒑₁∂γ(ap,ξ)

# ------------ Linear1D ---------------
@inline get𝑛𝒑(::ReproducingKernel{𝝃,:Linear1D}) where 𝝃 = 2
@inline get𝒑(::ReproducingKernel{𝝃,:Linear1D},x::NTuple{3,Float64}) where 𝝃 = (1.,x[1])
@inline get∂𝒑∂x(::ReproducingKernel{𝝃,:Linear1D},::NTuple{3,Float64}) where 𝝃 = (0.,1.)
@inline get∂𝒑∂y(::ReproducingKernel{𝝃,:Linear1D},::Any) where 𝝃 = (0.,0.)
@inline get∂𝒑∂z(::ReproducingKernel{𝝃,:Linear1D},::Any) where 𝝃 = (0.,0.)
@inline get∂²𝒑∂x²(::ReproducingKernel{𝝃,:Linear1D},::Any) where 𝝃 = (0.,0.)
@inline get∂²𝒑∂y²(::ReproducingKernel{𝝃,:Linear1D},::Any) where 𝝃 = (0.,0.)
@inline get∂²𝒑∂z²(::ReproducingKernel{𝝃,:Linear1D},::Any) where 𝝃 = (0.,0.)
@inline get∂²𝒑∂x∂y(::ReproducingKernel{𝝃,:Linear1D},::Any) where 𝝃 = (0.,0.)
@inline get∂²𝒑∂x∂z(::ReproducingKernel{𝝃,:Linear1D},::Any) where 𝝃 = (0.,0.)
@inline get∂²𝒑∂y∂z(::ReproducingKernel{𝝃,:Linear1D},::Any) where 𝝃 = (0.,0.)

@inline get𝑛𝒑₁(::ReproducingKernel{𝝃,:Linear1D}) where 𝝃 = 1
@inline get𝒑₁(::ReproducingKernel{𝝃,:Linear1D},::Any) where 𝝃<:AbstractNode = (1.0,)
@inline get∂𝒑₁∂ξ(::ReproducingKernel{𝝃,:Linear1D},::Any) where 𝝃<:AbstractNode = (0.0,)

# ------------ Quadaratic1D ---------------
@inline get𝑛𝒑(::ReproducingKernel{𝝃,:Quadratic1D}) where 𝝃 = 3
@inline get𝒑(::ReproducingKernel{𝝃,:Quadratic1D},x::NTuple{3,Float64}) where 𝝃 = (1.,x[1],x[1]^2)
@inline get∂𝒑∂x(::ReproducingKernel{𝝃,:Quadratic1D},x::NTuple{3,Float64}) where 𝝃 = (0.,1.,2*x[1])
@inline get∂𝒑∂y(::ReproducingKernel{𝝃,:Quadratic1D},::Any) where 𝝃 = (0.,0.,0.)
@inline get∂𝒑∂z(::ReproducingKernel{𝝃,:Quadratic1D},::Any) where 𝝃 = (0.,0.,0.)
@inline get∂²𝒑∂x²(::ReproducingKernel{𝝃,:Quadratic1D},::Any) where 𝝃 =(0.,0.,2.)
@inline get∂²𝒑∂y²(::ReproducingKernel{𝝃,:Quadratic1D},::Any) where 𝝃 =(0.,0.,0.)
@inline get∂²𝒑∂z²(::ReproducingKernel{𝝃,:Quadratic1D},::Any) where 𝝃 =(0.,0.,0.)
@inline get∂²𝒑∂x∂y(::ReproducingKernel{𝝃,:Quadratic1D},::Any) where 𝝃 =(0.,0.,0.)
@inline get∂²𝒑∂x∂z(::ReproducingKernel{𝝃,:Quadratic1D},::Any) where 𝝃 =(0.,0.,0.)
@inline get∂²𝒑∂y∂z(::ReproducingKernel{𝝃,:Quadratic1D},::Any) where 𝝃 =(0.,0.,0.)
@inline get∂³𝒑∂x³(::ReproducingKernel{𝝃,:Quadratic1D},::Any) where 𝝃 =(0.,0.,0.)
@inline get∂³𝒑∂x²∂y(::ReproducingKernel{𝝃,:Quadratic1D},::Any) where 𝝃 =(0.,0.,0.)
@inline get∂³𝒑∂x∂y²(::ReproducingKernel{𝝃,:Quadratic1D},::Any) where 𝝃 =(0.,0.,0.)
@inline get∂³𝒑∂y³(::ReproducingKernel{𝝃,:Quadratic1D},::Any) where 𝝃 =(0.,0.,0.)

@inline get𝑛𝒑₁(::ReproducingKernel{𝝃,:Quadratic1D}) where 𝝃 = 2
@inline get𝒑₁(ap::ReproducingKernel{𝝃,:Quadratic1D},ξ::𝝃) where 𝝃<:AbstractNode = get𝒑₁(ap,ξ.ξ)
@inline get𝒑₁(::ReproducingKernel{𝝃,:Quadratic1D},ξ::Float64) where 𝝃<:AbstractNode = (1.0,0.5*(1.0-ξ))
@inline get∂𝒑₁∂ξ(::ReproducingKernel{𝝃,:Quadratic1D},::Any) where 𝝃<:AbstractNode = (0.0,1.0)

# ------------ Cubic1D ---------------
@inline get𝑛𝒑(::ReproducingKernel{𝝃,:Cubic1D}) where 𝝃 = 4
@inline get𝒑(::ReproducingKernel{𝝃,:Cubic1D},x::NTuple{3,Float64}) where 𝝃 = (1.,x[1],x[1]^2,x[1]^3)
@inline get∂𝒑∂x(::ReproducingKernel{𝝃,:Cubic1D},x::NTuple{3,Float64}) where 𝝃 = (0.,1.,2*x[1],3*x[1]^2)
@inline get∂𝒑∂y(::ReproducingKernel{𝝃,:Cubic1D}, ::Any) where 𝝃 = (0.,0.,0.,0.)
@inline get∂𝒑∂z(::ReproducingKernel{𝝃,:Cubic1D}, ::Any) where 𝝃 = (0.,0.,0.,0.)
@inline get∂²𝒑∂x²(::ReproducingKernel{𝝃,:Cubic1D},x::NTuple{3,Float64}) where 𝝃 = (0.,0.,2.,6*x[1])
@inline get∂²𝒑∂y²(::ReproducingKernel{𝝃,:Cubic1D},x::NTuple{3,Float64}) where 𝝃 = (0.,0.,0.,0.)
@inline get∂²𝒑∂z²(::ReproducingKernel{𝝃,:Cubic1D},x::NTuple{3,Float64}) where 𝝃 = (0.,0.,0.,0.)
@inline get∂²𝒑∂x∂y(::ReproducingKernel{𝝃,:Cubic1D},x::NTuple{3,Float64}) where 𝝃 = (0.,0.,0.,0.)
@inline get∂²𝒑∂x∂z(::ReproducingKernel{𝝃,:Cubic1D},x::NTuple{3,Float64}) where 𝝃 = (0.,0.,0.,0.)
@inline get∂²𝒑∂y∂z(::ReproducingKernel{𝝃,:Cubic1D},x::NTuple{3,Float64}) where 𝝃 = (0.,0.,0.,0.)
@inline get∂³𝒑∂x³(::ReproducingKernel{𝝃,:Cubic1D},x::NTuple{3,Float64}) where 𝝃 = (0.,0.,0.,6.)
@inline get∂³𝒑∂x²∂y(::ReproducingKernel{𝝃,:Cubic1D},x::NTuple{3,Float64}) where 𝝃 = (0.,0.,0.,0.)
@inline get∂³𝒑∂x∂y²(::ReproducingKernel{𝝃,:Cubic1D},x::NTuple{3,Float64}) where 𝝃 = (0.,0.,0.,0.)
@inline get∂³𝒑∂y³(::ReproducingKernel{𝝃,:Cubic1D},x::NTuple{3,Float64}) where 𝝃 = (0.,0.,0.,0.)

@inline get𝑛𝒑₁(::ReproducingKernel{𝝃,:Cubic1D}) where 𝝃 = 3
@inline get𝒑₁(ap::ReproducingKernel{𝝃,:Cubic1D},ξ::𝝃) where 𝝃<:AbstractNode = get𝒑₁(ap,ξ.ξ)
@inline get𝒑₁(::ReproducingKernel{𝝃,:Cubic1D},ξ::Float64) where 𝝃<:AbstractNode = (1.0,0.5*(1.0-ξ),0.25*(1.0-ξ)^2)
@inline get∂𝒑₁∂ξ(ap::ReproducingKernel{𝝃,:Cubic1D},ξ::𝝃) where 𝝃<:AbstractNode = get∂𝒑₁∂ξ(ap,ξ.ξ)
@inline get∂𝒑₁∂ξ(::ReproducingKernel{𝝃,:Cubic1D},ξ::Float64) where 𝝃<:AbstractNode = (0.,1.0,(1.0-ξ))

# ------------ Linear2D ---------------
@inline get𝑛𝒑(::ReproducingKernel{𝝃,:Linear2D}) where 𝝃 = 3
@inline get𝒑(::ReproducingKernel{𝝃,:Linear2D},x::NTuple{3,Float64}) where 𝝃 = (1.,x[1],x[2])
@inline get∂𝒑∂x(::ReproducingKernel{𝝃,:Linear2D}, ::Any) where 𝝃 = (0.,1.,0.)
@inline get∂𝒑∂y(::ReproducingKernel{𝝃,:Linear2D}, ::Any) where 𝝃 = (0.,0.,1.)
@inline get∂𝒑∂z(::ReproducingKernel{𝝃,:Linear2D}, ::Any) where 𝝃 = (0.,0.,0.)

@inline get𝑛𝒑₁(::ReproducingKernel{𝝃,:Linear2D}) where 𝝃 = 1
@inline get𝒑₁(ap::ReproducingKernel{𝝃,:Linear2D},ξ::𝝃) where 𝝃<:AbstractNode = get𝒑₁(ap,ξ.ξ,ξ.η)
@inline get𝒑₁(::ReproducingKernel{𝝃,:Linear2D},::Any,::Any) where 𝝃<:AbstractNode = (1.,)
@inline get∂𝒑₁∂ξ(ap::ReproducingKernel{𝝃,:Linear2D},ξ::𝝃) where 𝝃<:AbstractNode = get∂𝒑₁∂ξ(ap,ξ.ξ,ξ.η)
@inline get∂𝒑₁∂ξ(::ReproducingKernel{𝝃,:Linear2D},::Any,::Any) where 𝝃<:AbstractNode = (0.,)
@inline get∂𝒑₁∂η(ap::ReproducingKernel{𝝃,:Linear2D},ξ::𝝃) where 𝝃<:AbstractNode = get∂𝒑₁∂η(ap,ξ.ξ,ξ.η)
@inline get∂𝒑₁∂η(::ReproducingKernel{𝝃,:Linear2D},::Any,::Any) where 𝝃<:AbstractNode = (0.,)

# ------------ Quadratic2D ---------------
@inline get𝑛𝒑(::ReproducingKernel{𝝃,:Quadratic2D}) where 𝝃 = 6
@inline get𝒑(::ReproducingKernel{𝝃,:Quadratic2D},x::NTuple{3,Float64}) where 𝝃 = (1.,x[1],x[2],x[1]^2,x[1]*x[2],x[2]^2)
@inline get∂𝒑∂x(::ReproducingKernel{𝝃,:Quadratic2D},x::NTuple{3,Float64}) where 𝝃 = (0.,1.,0.,2*x[1],x[2],0.)
@inline get∂𝒑∂y(::ReproducingKernel{𝝃,:Quadratic2D},x::NTuple{3,Float64}) where 𝝃 = (0.,0.,1.,0.,x[1],2*x[2])
@inline get∂𝒑∂z(::ReproducingKernel{𝝃,:Quadratic2D}, ::Any) where 𝝃 = (0.,0.,0.,0.,0.,0.)
@inline get∂²𝒑∂x²(::ReproducingKernel{𝝃,:Quadratic2D},x::NTuple{3,Float64}) where 𝝃 = (0.,0.,0.,2.,0.,0.)
@inline get∂²𝒑∂y²(::ReproducingKernel{𝝃,:Quadratic2D},x::NTuple{3,Float64}) where 𝝃 = (0.,0.,0.,0.,0.,2.)
@inline get∂²𝒑∂z²(::ReproducingKernel{𝝃,:Quadratic2D},x::NTuple{3,Float64}) where 𝝃 = (0.,0.,0.,0.,0.,0.)
@inline get∂²𝒑∂x∂y(::ReproducingKernel{𝝃,:Quadratic2D},x::NTuple{3,Float64}) where 𝝃 = (0.,0.,0.,0.,1.,0.)
@inline get∂²𝒑∂x∂z(::ReproducingKernel{𝝃,:Quadratic2D},x::NTuple{3,Float64}) where 𝝃 = (0.,0.,0.,0.,0.,0.)
@inline get∂²𝒑∂y∂z(::ReproducingKernel{𝝃,:Quadratic2D},x::NTuple{3,Float64}) where 𝝃 = (0.,0.,0.,0.,0.,0.)
@inline get∂³𝒑∂x³(::ReproducingKernel{𝝃,:Quadratic2D},x::NTuple{3,Float64}) where 𝝃 = (0.,0.,0.,0.,0.,0.)
@inline get∂³𝒑∂x²∂y(::ReproducingKernel{𝝃,:Quadratic2D},x::NTuple{3,Float64}) where 𝝃 = (0.,0.,0.,0.,0.,0.)
@inline get∂³𝒑∂x∂y²(::ReproducingKernel{𝝃,:Quadratic2D},x::NTuple{3,Float64}) where 𝝃 = (0.,0.,0.,0.,0.,0.)
@inline get∂³𝒑∂y³(::ReproducingKernel{𝝃,:Quadratic2D},x::NTuple{3,Float64}) where 𝝃 = (0.,0.,0.,0.,0.,0.)

@inline get𝑛𝒑₁(::ReproducingKernel{𝝃,:Quadratic2D}) where 𝝃 = 3
@inline get𝒑₁(ap::ReproducingKernel{𝝃,:Quadratic2D},ξ::𝝃) where 𝝃<:AbstractNode = get𝒑₁(ap,ξ.ξ,ξ.η)
@inline get𝒑₁(::ReproducingKernel{𝝃,:Quadratic2D},ξ::Float64,η::Float64) where 𝝃<:AbstractNode = (1.,ξ,η)
@inline get∂𝒑₁∂ξ(ap::ReproducingKernel{𝝃,:Quadratic2D},ξ::Any) where 𝝃 = (0.,1.,0.)
@inline get∂𝒑₁∂η(ap::ReproducingKernel{𝝃,:Quadratic2D},ξ::Any) where 𝝃 = (0.,0.,1.)
@inline get∂²𝒑₁∂ξ²(ap::ReproducingKernel{𝝃,:Quadratic2D},ξ::Any) where 𝝃 = (0.,0.,0.)
@inline get∂²𝒑₁∂ξ∂η(ap::ReproducingKernel{𝝃,:Quadratic2D},ξ::Any) where 𝝃 = (0.,0.,0.)
@inline get∂²𝒑₁∂η²(ap::ReproducingKernel{𝝃,:Quadratic2D},ξ::Any) where 𝝃 = (0.,0.,0.)
@inline get∂𝒑₁∂ξ(::ReproducingKernel{𝝃,:Quadratic2D},::Any,::Any) where 𝝃<:AbstractNode = (0.,1.,0.)
@inline get∂𝒑₁∂η(::ReproducingKernel{𝝃,:Quadratic2D},::Any,::Any) where 𝝃<:AbstractNode = (0.,0.,1.)

@inline get𝑛𝒑₂(::ReproducingKernel{𝝃,:Quadratic2D}) where 𝝃 = 1
@inline get𝒑₂(ap::ReproducingKernel{𝝃,:Quadratic2D},ξ::Any) where 𝝃 = (1.,)
@inline get∂𝒑₂∂ξ(ap::ReproducingKernel{𝝃,:Quadratic2D},ξ::Any) where 𝝃 = (0.,)
@inline get∂𝒑₂∂η(ap::ReproducingKernel{𝝃,:Quadratic2D},ξ::Any) where 𝝃 = (0.,)
@inline get∂²𝒑₂∂ξ²(ap::ReproducingKernel{𝝃,:Quadratic2D},ξ::Any) where 𝝃 = (0.,)
@inline get∂²𝒑₂∂ξ∂η(ap::ReproducingKernel{𝝃,:Quadratic2D},ξ::Any) where 𝝃 = (0.,)

# ------------ Cubic2D ---------------
@inline get𝑛𝒑(::ReproducingKernel{𝝃,:Cubic2D}) where 𝝃 = 10
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
@inline get∂²𝒑∂x²(::ReproducingKernel{𝝃,:Cubic2D},x::NTuple{3,Float64}) where 𝝃 =
(
    0., 0., 0., 2., 0., 0., 6*x[1], 2*x[2], 0., 0.
)
@inline get∂²𝒑∂x∂y(::ReproducingKernel{𝝃,:Cubic2D},x::NTuple{3,Float64}) where 𝝃 =
(
    0., 0., 0., 0., 1., 0., 0., 2*x[1], 2*x[2], 0.
)
@inline get∂²𝒑∂y²(::ReproducingKernel{𝝃,:Cubic2D},x::NTuple{3,Float64}) where 𝝃 =
(
    0., 0., 0., 0., 0., 2., 0., 0., 2*x[1], 6*x[2]
)
@inline get∂𝒑∂z(::ReproducingKernel{𝝃,:Cubic2D},::Any) where 𝝃 =
(
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
)
@inline get∂²𝒑∂x∂z(::ReproducingKernel{𝝃,:Cubic2D},x::NTuple{3,Float64}) where 𝝃 =
(
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
)
@inline get∂²𝒑∂y∂z(::ReproducingKernel{𝝃,:Cubic2D},x::NTuple{3,Float64}) where 𝝃 =
(
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
)
@inline get∂²𝒑∂z²(::ReproducingKernel{𝝃,:Cubic2D},x::NTuple{3,Float64}) where 𝝃 =
(
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
)
@inline get∂³𝒑∂x³(::ReproducingKernel{𝝃,:Cubic2D},x::NTuple{3,Float64}) where 𝝃 =
(
    0., 0., 0., 0., 0., 0., 6., 0., 0., 0.
)
@inline get∂³𝒑∂x²∂y(::ReproducingKernel{𝝃,:Cubic2D},x::NTuple{3,Float64}) where 𝝃 =
(
    0., 0., 0., 0., 0., 0., 0., 2., 0., 0.
)
@inline get∂³𝒑∂x∂y²(::ReproducingKernel{𝝃,:Cubic2D},x::NTuple{3,Float64}) where 𝝃 =
(
    0., 0., 0., 0., 0., 0., 0., 0., 2., 0.
)
@inline get∂³𝒑∂y³(::ReproducingKernel{𝝃,:Cubic2D},x::NTuple{3,Float64}) where 𝝃 =
(
    0., 0., 0., 0., 0., 0., 0., 0., 0., 6.
)

@inline get𝑛𝒑₁(::ReproducingKernel{𝝃,:Cubic2D}) where 𝝃 = 6
@inline get𝒑₁(ap::ReproducingKernel{𝝃,:Cubic2D},ξ::𝝃) where 𝝃<:AbstractNode = get𝒑₁(ap,ξ.ξ,ξ.η)
@inline get𝒑₁(::ReproducingKernel{𝝃,:Cubic2D},ξ::Float64,η::Float64) where 𝝃<:AbstractNode = (1.,ξ,η,ξ^2,ξ*η,η^2)
@inline get∂𝒑₁∂ξ(ap::ReproducingKernel{𝝃,:Cubic2D},ξ::𝝃) where 𝝃<:AbstractNode = get∂𝒑₁∂ξ(ap,ξ.ξ,ξ.η)
@inline get∂𝒑₁∂ξ(::ReproducingKernel{𝝃,:Cubic2D},ξ::Float64,η::Float64) where 𝝃<:AbstractNode = (0.,1.,0.,2.0*ξ,η,0.)
@inline get∂𝒑₁∂η(ap::ReproducingKernel{𝝃,:Cubic2D},ξ::𝝃) where 𝝃<:AbstractNode = get∂𝒑₁∂η(ap,ξ.ξ,ξ.η)
@inline get∂𝒑₁∂η(::ReproducingKernel{𝝃,:Cubic2D},ξ::Float64,η::Float64) where 𝝃<:AbstractNode = (0.,0.,1.,0.,ξ,2.0*η)

@inline get𝑛𝒑₂(::ReproducingKernel{𝝃,:Cubic2D}) where 𝝃 = 3
@inline get𝒑₂(ap::ReproducingKernel{𝝃,:Cubic2D},ξ::𝝃) where 𝝃<:AbstractNode = get𝒑₂(ap,ξ.ξ,ξ.η)
@inline get𝒑₂(::ReproducingKernel{𝝃,:Cubic2D},ξ::Float64,η::Float64) where 𝝃<:AbstractNode = (1.,ξ,η)
@inline get∂𝒑₂∂ξ(ap::ReproducingKernel{𝝃,:Cubic2D},ξ::𝝃) where 𝝃<:AbstractNode = (0.,1.,0.)
@inline get∂𝒑₂∂η(ap::ReproducingKernel{𝝃,:Cubic2D},ξ::𝝃) where 𝝃<:AbstractNode = (0.,0.,1.)
@inline get∂²𝒑₂∂ξ²(ap::ReproducingKernel{𝝃,:Cubic2D},ξ::Any) where 𝝃<:AbstractNode = (0.,0.,0.)
@inline get∂²𝒑₂∂ξ∂η(ap::ReproducingKernel{𝝃,:Cubic2D},ξ::Any) where 𝝃<:AbstractNode = (0.,0.,0.)
@inline get∂²𝒑₂∂η²(ap::ReproducingKernel{𝝃,:Cubic2D},ξ::Any) where 𝝃<:AbstractNode = (0.,0.,0.)

## Kernel Function
function get𝜙(ap::ReproducingKernel{𝝃,𝒑,:□,𝜙},x::Node,Δx::NTuple{3,Float64}) where {𝝃,𝒑,𝜙}
    rx = abs(Δx[1])/x.s₁
    ry = abs(Δx[2])/x.s₂
    rz = abs(Δx[3])/x.s₃
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

function get∇²𝜙(ap::ReproducingKernel{𝝃,𝒑,:□,𝜙},x::Node,Δx::NTuple{3,Float64}) where {𝝃,𝒑,𝜙}
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
    ∂²wx = get∂𝜙∂r(ap,rx)*∂rx^2
    ∂²wy = get∂𝜙∂r(ap,ry)*∂ry^2
    ∂²wz = get∂𝜙∂r(ap,rz)*∂rz^2
    return wx*wy*wz, ∂wx*wy*wz, wx*∂wy*wz, ∂²wx*wy*wz, ∂wx*∂wy*wz, wx*∂²wy*wz, wx*wy*∂wz, ∂wx*wy*∂wz, wx*∂wy*∂wz, wx*wy*∂²wz
end

function get∇³𝜙(ap::ReproducingKernel{𝝃,𝒑,:□,𝜙},x::Node,Δx::NTuple{3,Float64}) where {𝝃,𝒑,𝜙}
    rx = abs(Δx[1])/x.s₁
    ry = abs(Δx[2])/x.s₂
    ∂rx = sign(Δx[1])/x.s₁
    ∂ry = sign(Δx[2])/x.s₂
    wx = get𝜙ᵣ(ap,rx)
    wy = get𝜙ᵣ(ap,ry)
    ∂wx = get∂𝜙∂r(ap,rx)*∂rx
    ∂wy = get∂𝜙∂r(ap,ry)*∂ry
    ∂²wx = get∂𝜙∂r(ap,rx)*∂rx^2
    ∂²wy = get∂𝜙∂r(ap,ry)*∂ry^2
    ∂³wx = get∂𝜙∂r(ap,rx)*∂rx^3
    ∂³wy = get∂𝜙∂r(ap,ry)*∂ry^3
    return wx*wy, ∂wx*wy, wx*∂wy, ∂²wx*wy, ∂wx*∂wy, wx*∂²wy, ∂³wx*wy, ∂²wx*∂wy, ∂wx*∂²wy, wx*∂³wy
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

function get𝜙ᵣ(::ReproducingKernel{𝝃,𝒑,𝑠,:QuinticSpline},r::Float64) where {𝝃,𝒑,𝑠}
    if r > 1.0
        return 0.0
    elseif r <= 1/3
        return ((3-3r)^5 - 6(2-3r)^5 + 15(1-3r)^5)/120
    elseif r <= 2/3 && r > 1/3
        return ((3-3r)^5 - 6(2-3r)^5)/120
    else
        return (3-3r)^5/120
    end
end

function get∂𝜙∂r(::ReproducingKernel{𝝃,𝒑,𝑠,:QuinticSpline},r::Float64) where {𝝃,𝒑,𝑠}
    if r > 1.0
        return 0.0
    elseif r <= 1/3
        return -((3-3r)^4 - 6(2-3r)^4 + 15(1-3r)^4)/8
    elseif r <= 2/3 && r > 1/3
        return -((3-3r)^4 - 6(2-3r)^4)/8
    else
        return -(3-3r)^4/8
    end
end

function get∂²𝜙∂r²(::ReproducingKernel{𝝃,𝒑,𝑠,:QuinticSpline},r::Float64) where {𝝃,𝒑,𝑠}
    if r > 1.0
        return 0.0
    elseif r <= 1/3
        return ((3-3r)^3 - 6(2-3r)^3 + 15(1-3r)^3)*1.5
    elseif r <= 2/3 && r > 1/3
        return ((3-3r)^3 - 6(2-3r)^3)*1.5
    else
        return (3-3r)^3*1.5
    end
end

function get∂³𝜙∂r³(::ReproducingKernel{𝝃,𝒑,𝑠,:QuinticSpline},r::Float64) where {𝝃,𝒑,𝑠}
    if r > 1.0
        return 0.0
    elseif r <= 1/3
        return -((3-3r)^2 - 6(2-3r)^2 + 15(1-3r)^2)*13.5
    elseif r <= 2/3 && r > 1/3
        return -((3-3r)^2 - 6(2-3r)^2)*13.5
    else
        return -(3-3r)^2*13.5
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
        # print(Δx)
        𝒑 = get𝒑(ap,Δx)
        𝜙 = get𝜙(ap,xᵢ,Δx)
        # print(𝜙)
        for I in 1:n
            for J in I:n
                𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
            end
        end
    end
    # print(𝗠.m[1:55])
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

function cal∇²𝗠!(ap::ReproducingKernel,x::NTuple{3,Float64})
    𝓒 = ap.𝓒
    𝗠 = ap.𝗠[:∂1]
    ∂𝗠∂x = ap.𝗠[:∂x_]
    ∂𝗠∂y = ap.𝗠[:∂y_]
    ∂𝗠∂z = ap.𝗠[:∂z_]
    ∂𝗠⁻¹∂x = ap.𝗠[:∂x]
    ∂𝗠⁻¹∂y = ap.𝗠[:∂y]
    ∂𝗠⁻¹∂z = ap.𝗠[:∂z]
    ∂²𝗠∂x² = ap.𝗠[:∂x²]
    ∂²𝗠∂y² = ap.𝗠[:∂y²]
    ∂²𝗠∂z² = ap.𝗠[:∂z²]
    ∂²𝗠∂x∂y = ap.𝗠[:∂x∂y]
    ∂²𝗠∂x∂z = ap.𝗠[:∂x∂z]
    ∂²𝗠∂y∂z = ap.𝗠[:∂y∂z]
    n = get𝑛𝒑(ap)
    fill!(𝗠,0.)
    fill!(∂𝗠∂x,0.)
    fill!(∂𝗠∂y,0.)
    fill!(∂𝗠∂z,0.)
    fill!(∂²𝗠∂x²,0.)
    fill!(∂²𝗠∂y²,0.)
    fill!(∂²𝗠∂z²,0.)
    fill!(∂²𝗠∂x∂y,0.)
    fill!(∂²𝗠∂x∂z,0.)
    fill!(∂²𝗠∂y∂z,0.)
    for xᵢ in 𝓒
        Δx = x - xᵢ
        𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂²𝒑∂x², ∂²𝒑∂x∂y, ∂²𝒑∂y², ∂𝒑∂z, ∂²𝒑∂x∂z, ∂²𝒑∂y∂z, ∂²𝒑∂z² = get∇²𝒑(ap,Δx)
        𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂²𝜙∂x², ∂²𝜙∂x∂y, ∂²𝜙∂y², ∂𝜙∂z, ∂²𝜙∂x∂z, ∂²𝜙∂y∂z, ∂²𝜙∂z² = get∇²𝜙(ap,xᵢ,Δx)
        for I in 1:n
            for J in I:n
                𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
                ∂𝗠∂x[I,J] += ∂𝜙∂x*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂x[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂x[J]
                ∂𝗠∂y[I,J] += ∂𝜙∂y*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂y[J]
                ∂𝗠∂z[I,J] += ∂𝜙∂z*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂z[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂z[J]

                ∂²𝗠∂x²[I,J] += ∂²𝜙∂x²*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂x²[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂x²[J] + 2.0*∂𝜙∂x*∂𝒑∂x[I]*𝒑[J] + 2.0*∂𝜙∂x*𝒑[I]*∂𝒑∂x[J] + 2.0*𝜙*∂𝒑∂x[I]*∂𝒑∂x[J]

                ∂²𝗠∂y²[I,J] += ∂²𝜙∂y²*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂y²[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂y²[J] + 2.0*∂𝜙∂y*∂𝒑∂y[I]*𝒑[J] + 2.0*∂𝜙∂y*𝒑[I]*∂𝒑∂y[J] + 2.0*𝜙*∂𝒑∂y[I]*∂𝒑∂y[J]

                ∂²𝗠∂z²[I,J] += ∂²𝜙∂z²*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂z²[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂z²[J] + 2.0*∂𝜙∂z*∂𝒑∂z[I]*𝒑[J] + 2.0*∂𝜙∂z*𝒑[I]*∂𝒑∂z[J] + 2.0*𝜙*∂𝒑∂z[I]*∂𝒑∂z[J]

                ∂²𝗠∂x∂y[I,J] += ∂²𝜙∂x∂y*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂x∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂x∂y[J] + ∂𝜙∂x*∂𝒑∂y[I]*𝒑[J] + ∂𝜙∂y*∂𝒑∂x[I]*𝒑[J] + ∂𝜙∂x*𝒑[I]*∂𝒑∂y[J] + ∂𝜙∂y*𝒑[I]*∂𝒑∂x[J] + 𝜙*∂𝒑∂x[I]*∂𝒑∂y[J] + 𝜙*∂𝒑∂y[I]*∂𝒑∂x[J]

                ∂²𝗠∂x∂z[I,J] += ∂²𝜙∂x∂z*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂x∂z[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂x∂z[J] + ∂𝜙∂x*∂𝒑∂z[I]*𝒑[J] + ∂𝜙∂z*∂𝒑∂x[I]*𝒑[J] + ∂𝜙∂x*𝒑[I]*∂𝒑∂z[J] + ∂𝜙∂z*𝒑[I]*∂𝒑∂x[J] + 𝜙*∂𝒑∂x[I]*∂𝒑∂z[J] + 𝜙*∂𝒑∂z[I]*∂𝒑∂x[J]

                ∂²𝗠∂y∂z[I,J] += ∂²𝜙∂y∂z*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂y∂z[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂y∂z[J] + ∂𝜙∂y*∂𝒑∂z[I]*𝒑[J] + ∂𝜙∂z*∂𝒑∂y[I]*𝒑[J] + ∂𝜙∂y*𝒑[I]*∂𝒑∂z[J] + ∂𝜙∂z*𝒑[I]*∂𝒑∂y[J] + 𝜙*∂𝒑∂y[I]*∂𝒑∂z[J] + 𝜙*∂𝒑∂z[I]*∂𝒑∂y[J]
            end
        end
    end
    cholesky!(𝗠)
    U⁻¹ = inverse!(𝗠)
    ∂𝗠⁻¹∂x = - UUᵀAUUᵀ!(∂𝗠⁻¹∂x,∂𝗠∂x,U⁻¹)
    ∂𝗠⁻¹∂y = - UUᵀAUUᵀ!(∂𝗠⁻¹∂y,∂𝗠∂y,U⁻¹)
    ∂𝗠⁻¹∂z = - UUᵀAUUᵀ!(∂𝗠⁻¹∂z,∂𝗠∂z,U⁻¹)
    ∂²𝗠⁻¹∂x² = UUᵀAUUᵀ!(∂²𝗠∂x²,U⁻¹)
    ∂²𝗠⁻¹∂y² = UUᵀAUUᵀ!(∂²𝗠∂y²,U⁻¹)
    ∂²𝗠⁻¹∂z² = UUᵀAUUᵀ!(∂²𝗠∂z²,U⁻¹)
    ∂²𝗠⁻¹∂x∂y = UUᵀAUUᵀ!(∂²𝗠∂x∂y,U⁻¹)
    ∂²𝗠⁻¹∂x∂z = UUᵀAUUᵀ!(∂²𝗠∂x∂z,U⁻¹)
    ∂²𝗠⁻¹∂y∂z = UUᵀAUUᵀ!(∂²𝗠∂y∂z,U⁻¹)
    𝗠⁻¹ = UUᵀ!(U⁻¹)
    for i in 1:n
        for j in i:n
            for k in 1:n
                for l in 1:n
                    ∂²𝗠⁻¹∂x²[i,j] += 2*𝗠⁻¹[i,k]*∂𝗠∂x[k,l]*∂𝗠⁻¹∂x[l,j]
                    ∂²𝗠⁻¹∂y²[i,j] += 2*𝗠⁻¹[i,k]*∂𝗠∂y[k,l]*∂𝗠⁻¹∂y[l,j]
                    ∂²𝗠⁻¹∂z²[i,j] += 2*𝗠⁻¹[i,k]*∂𝗠∂z[k,l]*∂𝗠⁻¹∂z[l,j]
                    ∂²𝗠⁻¹∂x∂y[i,j] += 𝗠⁻¹[i,k]*(∂𝗠∂x[k,l]*∂𝗠⁻¹∂y[l,j] + ∂𝗠∂y[k,l]*∂𝗠⁻¹∂x[l,j])
                    ∂²𝗠⁻¹∂x∂z[i,j] += 𝗠⁻¹[i,k]*(∂𝗠∂x[k,l]*∂𝗠⁻¹∂z[l,j] + ∂𝗠∂z[k,l]*∂𝗠⁻¹∂x[l,j])
                    ∂²𝗠⁻¹∂y∂z[i,j] += 𝗠⁻¹[i,k]*(∂𝗠∂y[k,l]*∂𝗠⁻¹∂z[l,j] + ∂𝗠∂z[k,l]*∂𝗠⁻¹∂y[l,j])
                end
            end
        end
    end
    ∂²𝗠⁻¹∂x² = - ∂²𝗠⁻¹∂x²
    ∂²𝗠⁻¹∂y² = - ∂²𝗠⁻¹∂y²
    ∂²𝗠⁻¹∂z² = - ∂²𝗠⁻¹∂z²
    ∂²𝗠⁻¹∂x∂y = - ∂²𝗠⁻¹∂x∂y
    ∂²𝗠⁻¹∂x∂z = - ∂²𝗠⁻¹∂x∂z
    ∂²𝗠⁻¹∂y∂z = - ∂²𝗠⁻¹∂y∂z
    return 𝗠⁻¹, ∂𝗠⁻¹∂x, ∂𝗠⁻¹∂y, ∂²𝗠⁻¹∂x², ∂²𝗠⁻¹∂x∂y, ∂²𝗠⁻¹∂y², ∂𝗠⁻¹∂z, ∂²𝗠⁻¹∂x∂z, ∂²𝗠⁻¹∂y∂z, ∂²𝗠⁻¹∂z²
end

function cal∇³𝗠!(ap::ReproducingKernel,x::NTuple{3,Float64})
    𝓒 = ap.𝓒
    𝗠 = ap.𝗠[:∂1]
    ∂𝗠∂x = ap.𝗠[:∂x_]
    ∂𝗠∂y = ap.𝗠[:∂y_]
    ∂𝗠⁻¹∂x = ap.𝗠[:∂x]
    ∂𝗠⁻¹∂y = ap.𝗠[:∂y]
    ∂²𝗠∂x² = ap.𝗠[:∂x²_]
    ∂²𝗠∂x∂y = ap.𝗠[:∂x∂y_]
    ∂²𝗠∂y² = ap.𝗠[:∂y²_]
    ∂²𝗠⁻¹∂x² = ap.𝗠[:∂x²]
    ∂²𝗠⁻¹∂x∂y = ap.𝗠[:∂x∂y]
    ∂²𝗠⁻¹∂y² = ap.𝗠[:∂y²]
    ∂³𝗠∂x³ = ap.𝗠[:∂x³]
    ∂³𝗠∂x²∂y = ap.𝗠[:∂x²∂y]
    ∂³𝗠∂x∂y² = ap.𝗠[:∂x∂y²]
    ∂³𝗠∂y³ = ap.𝗠[:∂y³]
    n = get𝑛𝒑(ap)
    fill!(𝗠,0.)
    fill!(∂𝗠∂x,0.)
    fill!(∂𝗠∂y,0.)
    fill!(∂²𝗠∂x²,0.)
    fill!(∂²𝗠∂x∂y,0.)
    fill!(∂²𝗠∂y²,0.)
    fill!(∂³𝗠∂x³,0.)
    fill!(∂³𝗠∂x²∂y,0.)
    fill!(∂³𝗠∂x∂y²,0.)
    fill!(∂³𝗠∂y³,0.)
    for xᵢ in 𝓒
        Δx = x - xᵢ
        𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂²𝒑∂x², ∂²𝒑∂x∂y, ∂²𝒑∂y², ∂³𝒑∂x³, ∂³𝒑∂x²∂y, ∂³𝒑∂x∂y², ∂³𝒑∂y³ = get∇³𝒑(ap,Δx)
        𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂²𝜙∂x², ∂²𝜙∂x∂y, ∂²𝜙∂y², ∂³𝜙∂x³, ∂³𝜙∂x²∂y, ∂³𝜙∂x∂y², ∂³𝜙∂y³ = get∇³𝜙(ap,xᵢ,Δx)
        for I in 1:n
            for J in I:n
                𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
                ∂𝗠∂x[I,J] += ∂𝜙∂x*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂x[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂x[J]
                ∂𝗠∂y[I,J] += ∂𝜙∂y*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂y[J]

                ∂²𝗠∂x²[I,J] += ∂²𝜙∂x²*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂x²[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂x²[J] + 2.0*∂𝜙∂x*∂𝒑∂x[I]*𝒑[J] + 2.0*∂𝜙∂x*𝒑[I]*∂𝒑∂x[J] + 2.0*𝜙*∂𝒑∂x[I]*∂𝒑∂x[J]

                ∂²𝗠∂x∂y[I,J] += ∂²𝜙∂x∂y*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂x∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂x∂y[J] + ∂𝜙∂x*∂𝒑∂y[I]*𝒑[J] + ∂𝜙∂y*∂𝒑∂x[I]*𝒑[J] + ∂𝜙∂x*𝒑[I]*∂𝒑∂y[J] + ∂𝜙∂y*𝒑[I]*∂𝒑∂x[J] + 𝜙*∂𝒑∂x[I]*∂𝒑∂y[J] + 𝜙*∂𝒑∂y[I]*∂𝒑∂x[J]

                ∂²𝗠∂y²[I,J] += ∂²𝜙∂y²*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂y²[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂y²[J] + 2.0*∂𝜙∂y*∂𝒑∂y[I]*𝒑[J] + 2.0*∂𝜙∂y*𝒑[I]*∂𝒑∂y[J] + 2.0*𝜙*∂𝒑∂y[I]*∂𝒑∂y[J]

                ∂³𝗠∂x³[I,J] += ∂³𝜙∂x³*𝒑[I]*𝒑[J] + 𝜙*∂³𝒑∂x³[I]*𝒑[J] + 𝜙*𝒑[I]*∂³𝒑∂x³[J] + 3.0*∂²𝜙∂x²*∂𝒑∂x[I]*𝒑[J] + 3.0*∂𝜙∂x*∂²𝒑∂x²[I]*𝒑[J] + 3.0*∂²𝜙∂x²*𝒑[I]*∂𝒑∂x[J] + 3.0*∂𝜙∂x*𝒑[I]*∂²𝒑∂x²[J] + 3.0*𝜙*∂²𝒑∂x²[I]*∂𝒑∂x[J] + 3.0*𝜙*∂𝒑∂x[I]*∂²𝒑∂x²[J] + 6.0*∂𝜙∂x*∂𝒑∂x[I]*∂𝒑∂x[J]

                ∂³𝗠∂x²∂y[I,J] += ∂³𝜙∂x²∂y*𝒑[I]*𝒑[J] + 𝜙*∂³𝒑∂x²∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂³𝒑∂x²∂y[J] + 2.0*∂²𝜙∂x∂y*∂𝒑∂x[I]*𝒑[J] + ∂²𝜙∂x²*∂𝒑∂y[I]*𝒑[J] + 2.0*∂𝜙∂x*∂²𝒑∂x∂y[I]*𝒑[J] + ∂𝜙∂y*∂²𝒑∂x²[I]*𝒑[J] + 2.0*∂²𝜙∂x∂y*𝒑[I]*∂𝒑∂x[J] + ∂²𝜙∂x²*𝒑[I]*∂𝒑∂y[J] + 2.0*∂𝜙∂x*𝒑[I]*∂²𝒑∂x∂y[J] + ∂𝜙∂y*𝒑[I]*∂²𝒑∂x²[J] + 2.0*𝜙*∂²𝒑∂x∂y[I]*∂𝒑∂x[J] + 𝜙*∂²𝒑∂x²[I]*∂𝒑∂y[J] + 2.0*𝜙*∂𝒑∂x[I]*∂²𝒑∂x∂y[J] + 𝜙*∂𝒑∂y[I]*∂²𝒑∂x²[J] + 2.0*∂𝜙∂y*∂𝒑∂x[I]*∂𝒑∂x[J] + 2.0*∂𝜙∂x*∂𝒑∂y[I]*∂𝒑∂x[J] + 2.0*∂𝜙∂x*∂𝒑∂x[I]*∂𝒑∂y[J]

                ∂³𝗠∂x∂y²[I,J] += ∂³𝜙∂x∂y²*𝒑[I]*𝒑[J] + 𝜙*∂³𝒑∂x∂y²[I]*𝒑[J] + 𝜙*𝒑[I]*∂³𝒑∂x∂y²[J] + 2.0*∂²𝜙∂x∂y*∂𝒑∂y[I]*𝒑[J] + ∂²𝜙∂y²*∂𝒑∂x[I]*𝒑[J] + 2.0*∂𝜙∂y*∂²𝒑∂x∂y[I]*𝒑[J] + ∂𝜙∂x*∂²𝒑∂y²[I]*𝒑[J] + 2.0*∂²𝜙∂x∂y*𝒑[I]*∂𝒑∂y[J] + ∂²𝜙∂y²*𝒑[I]*∂𝒑∂x[J] + 2.0*∂𝜙∂y*𝒑[I]*∂²𝒑∂x∂y[J] + ∂𝜙∂x*𝒑[I]*∂²𝒑∂y²[J] + 2.0*𝜙*∂²𝒑∂x∂y[I]*∂𝒑∂y[J] + 𝜙*∂²𝒑∂y²[I]*∂𝒑∂x[J] + 2.0*𝜙*∂𝒑∂y[I]*∂²𝒑∂x∂y[J] + 𝜙*∂𝒑∂x[I]*∂²𝒑∂y²[J] + 2.0*∂𝜙∂x*∂𝒑∂y[I]*∂𝒑∂y[J] + 2.0*∂𝜙∂y*∂𝒑∂x[I]*∂𝒑∂y[J] + 2.0*∂𝜙∂y*∂𝒑∂y[I]*∂𝒑∂x[J]

                ∂³𝗠∂y³[I,J] += ∂³𝜙∂y³*𝒑[I]*𝒑[J] + 𝜙*∂³𝒑∂y³[I]*𝒑[J] + 𝜙*𝒑[I]*∂³𝒑∂y³[J] + 3.0*∂²𝜙∂y²*∂𝒑∂y[I]*𝒑[J] + 3.0*∂𝜙∂y*∂²𝒑∂y²[I]*𝒑[J] + 3.0*∂²𝜙∂y²*𝒑[I]*∂𝒑∂y[J] + 3.0*∂𝜙∂y*𝒑[I]*∂²𝒑∂y²[J] + 3.0*𝜙*∂²𝒑∂y²[I]*∂𝒑∂y[J] + 3.0*𝜙*∂𝒑∂y[I]*∂²𝒑∂y²[J] + 6.0*∂𝜙∂y*∂𝒑∂y[I]*∂𝒑∂y[J]
            end
        end
    end
    cholesky!(𝗠)
    U⁻¹ = inverse!(𝗠)
    ∂𝗠⁻¹∂x = - UUᵀAUUᵀ!(∂𝗠⁻¹∂x,∂𝗠∂x,U⁻¹)
    ∂𝗠⁻¹∂y = - UUᵀAUUᵀ!(∂𝗠⁻¹∂y,∂𝗠∂y,U⁻¹)
    ∂²𝗠⁻¹∂x² = UUᵀAUUᵀ!(∂²𝗠⁻¹∂x²,∂²𝗠∂x²,U⁻¹)
    ∂²𝗠⁻¹∂x∂y = UUᵀAUUᵀ!(∂²𝗠⁻¹∂x∂y,∂²𝗠∂x∂y,U⁻¹)
    ∂²𝗠⁻¹∂y² = UUᵀAUUᵀ!(∂²𝗠⁻¹∂y²,∂²𝗠∂y²,U⁻¹)
    ∂³𝗠⁻¹∂x³ = UUᵀAUUᵀ!(∂³𝗠∂x³,U⁻¹)
    ∂³𝗠⁻¹∂x²∂y = UUᵀAUUᵀ!(∂³𝗠∂x²∂y,U⁻¹)
    ∂³𝗠⁻¹∂x∂y² = UUᵀAUUᵀ!(∂³𝗠∂x∂y²,U⁻¹)
    ∂³𝗠⁻¹∂y³ = UUᵀAUUᵀ!(∂³𝗠∂y³,U⁻¹)
    𝗠⁻¹ = UUᵀ!(U⁻¹)
    for i in 1:n
        for j in i:n
            for k in 1:n
                for l in 1:n
                    ∂²𝗠⁻¹∂x²[i,j] += 2*𝗠⁻¹[i,k]*∂𝗠∂x[k,l]*∂𝗠⁻¹∂x[l,j]
                    ∂²𝗠⁻¹∂x∂y[i,j] += 𝗠⁻¹[i,k]*(∂𝗠∂x[k,l]*∂𝗠⁻¹∂y[l,j] + ∂𝗠∂y[k,l]*∂𝗠⁻¹∂x[l,j])
                    ∂²𝗠⁻¹∂y²[i,j] += 2*𝗠⁻¹[i,k]*∂𝗠∂y[k,l]*∂𝗠⁻¹∂y[l,j]
                end
            end
        end
    end
    ∂²𝗠⁻¹∂x² = - ∂²𝗠⁻¹∂x²
    ∂²𝗠⁻¹∂x∂y = - ∂²𝗠⁻¹∂x∂y
    ∂²𝗠⁻¹∂y² = - ∂²𝗠⁻¹∂y²
    for i in 1:n
        for j in i:n
            for k in 1:n
                for l in 1:n
                    ∂³𝗠⁻¹∂x³[i,j] += 3*𝗠⁻¹[i,k]*(∂²𝗠∂x²[k,l]*∂𝗠⁻¹∂x[l,j] + ∂𝗠∂x[k,l]*∂²𝗠⁻¹∂x²[l,j])
                    ∂³𝗠⁻¹∂x²∂y[i,j] += 𝗠⁻¹[i,k]*(2*∂²𝗠∂x∂y[k,l]*∂𝗠⁻¹∂x[l,j] + ∂²𝗠∂x²[k,l]*∂𝗠⁻¹∂y[l,j] + 2*∂𝗠∂x[k,l]*∂²𝗠⁻¹∂x∂y[l,j] + ∂𝗠∂y[k,l]*∂²𝗠⁻¹∂x²[l,j])
                    ∂³𝗠⁻¹∂x∂y²[i,j] += 𝗠⁻¹[i,k]*(2*∂²𝗠∂x∂y[k,l]*∂𝗠⁻¹∂y[l,j] + ∂²𝗠∂y²[k,l]*∂𝗠⁻¹∂x[l,j] + 2*∂𝗠∂y[k,l]*∂²𝗠⁻¹∂x∂y[l,j] + ∂𝗠∂x[k,l]*∂²𝗠⁻¹∂y²[l,j])
                    ∂³𝗠⁻¹∂y³[i,j] += 3*𝗠⁻¹[i,k]*(∂²𝗠∂y²[k,l]*∂𝗠⁻¹∂y[l,j] + ∂𝗠∂y[k,l]*∂²𝗠⁻¹∂y²[l,j])
                end
            end
        end
    end
    ∂³𝗠⁻¹∂x³ = - ∂³𝗠⁻¹∂x³
    ∂³𝗠⁻¹∂x²∂y = - ∂³𝗠⁻¹∂x²∂y
    ∂³𝗠⁻¹∂x∂y² = - ∂³𝗠⁻¹∂x∂y²
    ∂³𝗠⁻¹∂y³ = - ∂³𝗠⁻¹∂y³
    return 𝗠⁻¹, ∂𝗠⁻¹∂x, ∂𝗠⁻¹∂y, ∂²𝗠⁻¹∂x², ∂²𝗠⁻¹∂x∂y, ∂²𝗠⁻¹∂y², ∂³𝗠⁻¹∂x³, ∂³𝗠⁻¹∂x²∂y, ∂³𝗠⁻¹∂x∂y², ∂³𝗠⁻¹∂y³
end

function cal𝗚!(ap::ReproducingKernel)
    𝓖 = ap.𝓖
    𝗚 = ap.𝗠[:∇̃]
    n = get𝑛𝒒(ap)
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

function cal𝗚!(ap::ReproducingKernel{𝝃,:Quadratic1D,𝑠,𝜙,:Seg2}) where {𝝃<:AbstractNode,𝑠,𝜙}
    𝗚⁻¹ = ap.𝗠[:∇̃]
    fill!(𝗚⁻¹,0.0)
    𝐿 = get𝐿(ap)
    𝗚⁻¹[1] =  4.0/𝐿
    𝗚⁻¹[2] = -6.0/𝐿
    𝗚⁻¹[3] = 12.0/𝐿
    return 𝗚⁻¹
end

function cal𝗚!(ap::ReproducingKernel{𝝃,:Cubic1D,𝑠,𝜙,:Seg2}) where {𝝃<:AbstractNode,𝑠,𝜙}
    𝗚⁻¹ = ap.𝗠[:∇̃]
    fill!(𝗚⁻¹,0.0)
    𝐿 = get𝐿(ap)
    𝗚⁻¹[1] =    9.0/𝐿
    𝗚⁻¹[2] =  -36.0/𝐿
    𝗚⁻¹[3] =  192.0/𝐿
    𝗚⁻¹[4] =   30.0/𝐿
    𝗚⁻¹[5] = -180.0/𝐿
    𝗚⁻¹[6] =  180.0/𝐿
    return 𝗚⁻¹
end

function cal𝗚!(ap::ReproducingKernel{𝝃,:Linear2D,𝑠,𝜙,:Tri3}) where {𝝃<:AbstractNode,𝑠,𝜙}
    𝗚⁻¹ = ap.𝗠[:∇̃]
    fill!(𝗚⁻¹,0.0)
    𝐴 = get𝐴(ap)
    𝗚⁻¹[1] = 1.0/𝐴
    return 𝗚⁻¹
end

function cal𝗚!(ap::ReproducingKernel{𝝃,:Quadratic2D,𝑠,𝜙,:Tri3}) where {𝝃<:AbstractNode,𝑠,𝜙}
    𝗚⁻¹ = ap.𝗠[:∇̃]
    fill!(𝗚⁻¹,0.0)
    𝐴 = get𝐴(ap)
    𝗚⁻¹[1] =   9.0/𝐴
    𝗚⁻¹[2] = -12.0/𝐴
    𝗚⁻¹[3] =  24.0/𝐴
    𝗚⁻¹[4] = -12.0/𝐴
    𝗚⁻¹[5] =  12.0/𝐴
    𝗚⁻¹[6] =  24.0/𝐴
    return 𝗚⁻¹
end

function cal𝗚!(ap::ReproducingKernel{𝝃,:Cubic2D,𝑠,𝜙,:Tri3}) where {𝝃<:AbstractNode,𝑠,𝜙}
    𝗚⁻¹ = ap.𝗠[:∇̃]
    fill!(𝗚⁻¹,0.0)
    𝐴 = get𝐴(ap)
    𝗚⁻¹[1] =   36.0/𝐴
    𝗚⁻¹[2] = -120.0/𝐴
    𝗚⁻¹[3] =  600.0/𝐴
    𝗚⁻¹[4] = -120.0/𝐴
    𝗚⁻¹[5] =  300.0/𝐴
    𝗚⁻¹[6] =  600.0/𝐴
    𝗚⁻¹[7] =   90.0/𝐴
    𝗚⁻¹[8] = -540.0/𝐴
    𝗚⁻¹[9] = -180.0/𝐴
    𝗚⁻¹[10] =  540.0/𝐴
    𝗚⁻¹[11] =  180.0/𝐴
    𝗚⁻¹[12] = -720.0/𝐴
    𝗚⁻¹[13] = -720.0/𝐴
    𝗚⁻¹[14] =  540.0/𝐴
    𝗚⁻¹[15] = 1440.0/𝐴
    𝗚⁻¹[16] =   90.0/𝐴
    𝗚⁻¹[17] = -180.0/𝐴
    𝗚⁻¹[18] = -540.0/𝐴
    𝗚⁻¹[19] =   90.0/𝐴
    𝗚⁻¹[20] =  540.0/𝐴
    𝗚⁻¹[21] =  540.0/𝐴
    return 𝗚⁻¹
end

function cal𝗚₂!(ap::ReproducingKernel{𝝃,:Quadratic2D,𝑠,𝜙,:Tri3}) where {𝝃<:AbstractNode,𝑠,𝜙}
    𝗚⁻¹ = ap.𝗠[:∇̃]
    fill!(𝗚⁻¹,0.0)
    𝐴 = get𝐴(ap)
    𝗚⁻¹[1] = 1.0/𝐴
    return 𝗚⁻¹
end
function cal𝗚₂!(ap::ReproducingKernel{𝝃,:Cubic2D,𝑠,𝜙,:Tri3}) where {𝝃<:AbstractNode,𝑠,𝜙}
    𝗚⁻¹ = ap.𝗠[:∇̃]
    fill!(𝗚⁻¹,0.0)
    𝐴 = get𝐴(ap)
    𝗚⁻¹[1] =   9.0/𝐴
    𝗚⁻¹[2] = -12.0/𝐴
    𝗚⁻¹[3] =  24.0/𝐴
    𝗚⁻¹[4] = -12.0/𝐴
    𝗚⁻¹[5] =  12.0/𝐴
    𝗚⁻¹[6] =  24.0/𝐴
    return 𝗚⁻¹
end

function cal𝗚₂!(ap::ReproducingKernel{𝝃,:Quartic2D,𝑠,𝜙,:Tri3}) where {𝝃<:AbstractNode,𝑠,𝜙}
    𝗚⁻¹ = ap.𝗠[:∇̃]
    fill!(𝗚⁻¹,0.0)
    𝐴 = get𝐴(ap)
    𝗚⁻¹[1] =   36.0/𝐴
    𝗚⁻¹[2] = -120.0/𝐴
    𝗚⁻¹[3] =  600.0/𝐴
    𝗚⁻¹[4] = -120.0/𝐴
    𝗚⁻¹[5] =  300.0/𝐴
    𝗚⁻¹[6] =  600.0/𝐴
    𝗚⁻¹[7] =   90.0/𝐴
    𝗚⁻¹[8] = -540.0/𝐴
    𝗚⁻¹[9] = -180.0/𝐴
    𝗚⁻¹[10] =  540.0/𝐴
    𝗚⁻¹[11] =  180.0/𝐴
    𝗚⁻¹[12] = -720.0/𝐴
    𝗚⁻¹[13] = -720.0/𝐴
    𝗚⁻¹[14] =  540.0/𝐴
    𝗚⁻¹[15] = 1440.0/𝐴
    𝗚⁻¹[16] =   90.0/𝐴
    𝗚⁻¹[17] = -180.0/𝐴
    𝗚⁻¹[18] = -540.0/𝐴
    𝗚⁻¹[19] =   90.0/𝐴
    𝗚⁻¹[20] =  540.0/𝐴
    𝗚⁻¹[21] =  540.0/𝐴
    return 𝗚⁻¹
end
## shape functions
function get𝝭(ap::ReproducingKernel,ξ::Node)
    𝒙 = get𝒙(ap,ξ)
    return get𝝭(ap,𝒙)
end
function get∇𝝭(ap::ReproducingKernel,ξ::Node)
    𝒙 = get𝒙(ap,ξ)
    return get∇𝝭(ap,𝒙)
end
function get∇²𝝭(ap::ReproducingKernel,ξ::Node)
    𝒙 = get𝒙(ap,ξ)
    return get∇²𝝭(ap,𝒙)
end
function get∇³𝝭(ap::ReproducingKernel,ξ::Node)
    𝒙 = get𝒙(ap,ξ)
    return get∇³𝝭(ap,𝒙)
end

function get𝝭(ap::ReproducingKernel,𝒙::NTuple{3,Float64})
    𝓒 = ap.𝓒
    𝝭 = ap.𝝭[:∂1]
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

function get∇𝝭(ap::ReproducingKernel,𝒙::NTuple{3,Float64})
    𝓒 = ap.𝓒
    𝝭 = ap.𝝭[:∂1]
    ∂𝝭∂x = ap.𝝭[:∂x]
    ∂𝝭∂y = ap.𝝭[:∂y]
    ∂𝝭∂z = ap.𝝭[:∂z]
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

function get∇²𝝭(ap::ReproducingKernel,𝒙::NTuple{3,Float64})
    𝓒 = ap.𝓒
    𝝭 = ap.𝝭[:∂1]
    ∂𝝭∂x = ap.𝝭[:∂x]
    ∂𝝭∂y = ap.𝝭[:∂y]
    ∂𝝭∂z = ap.𝝭[:∂z]
    ∂²𝝭∂x² = ap.𝝭[:∂x²]
    ∂²𝝭∂y² = ap.𝝭[:∂y²]
    ∂²𝝭∂z² = ap.𝝭[:∂z²]
    ∂²𝝭∂x∂y = ap.𝝭[:∂x∂y]
    ∂²𝝭∂x∂z = ap.𝝭[:∂x∂z]
    ∂²𝝭∂y∂z = ap.𝝭[:∂y∂z]
    𝒑₀ᵀ𝗠⁻¹, 𝒑₀ᵀ∂𝗠⁻¹∂x, 𝒑₀ᵀ∂𝗠⁻¹∂y, 𝒑₀ᵀ∂²𝗠⁻¹∂x², 𝒑₀ᵀ∂²𝗠⁻¹∂x∂y, 𝒑₀ᵀ∂²𝗠⁻¹∂y², 𝒑₀ᵀ∂𝗠⁻¹∂z, 𝒑₀ᵀ∂²𝗠⁻¹∂x∂z, 𝒑₀ᵀ∂²𝗠⁻¹∂y∂z, 𝒑₀ᵀ∂²𝗠⁻¹∂z² = cal∇²𝗠!(ap,𝒙)
    for i in 1:length(𝓒)
        𝒙ᵢ = 𝓒[i]
        Δ𝒙 = 𝒙 - 𝒙ᵢ
        𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂²𝒑∂x², ∂²𝒑∂x∂y, ∂²𝒑∂y², ∂𝒑∂z, ∂²𝒑∂x∂z, ∂²𝒑∂y∂z, ∂²𝒑∂z² = get∇²𝒑(ap,Δ𝒙)
        𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂²𝜙∂x², ∂²𝜙∂x∂y, ∂²𝜙∂y², ∂𝜙∂z, ∂²𝜙∂x∂z, ∂²𝜙∂y∂z, ∂²𝜙∂z² = get∇²𝜙(ap,𝒙ᵢ,Δ𝒙)
        𝒑₀ᵀ𝗠⁻¹𝒑 = 𝒑₀ᵀ𝗠⁻¹*𝒑
        𝒑₀ᵀ∂𝗠⁻¹∂x𝒑 = 𝒑₀ᵀ∂𝗠⁻¹∂x*𝒑
        𝒑₀ᵀ∂𝗠⁻¹∂y𝒑 = 𝒑₀ᵀ∂𝗠⁻¹∂y*𝒑
        𝒑₀ᵀ∂𝗠⁻¹∂z𝒑 = 𝒑₀ᵀ∂𝗠⁻¹∂z*𝒑
        𝒑₀ᵀ𝗠⁻¹∂𝒑∂x = 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂x
        𝒑₀ᵀ𝗠⁻¹∂𝒑∂y = 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂y
        𝒑₀ᵀ𝗠⁻¹∂𝒑∂z = 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂z
        𝒑₀ᵀ∂²𝗠⁻¹∂x²𝒑 = 𝒑₀ᵀ∂²𝗠⁻¹∂x²*𝒑
        𝒑₀ᵀ∂²𝗠⁻¹∂y²𝒑 = 𝒑₀ᵀ∂²𝗠⁻¹∂y²*𝒑
        𝒑₀ᵀ∂²𝗠⁻¹∂z²𝒑 = 𝒑₀ᵀ∂²𝗠⁻¹∂z²*𝒑
        𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x² = 𝒑₀ᵀ𝗠⁻¹*∂²𝒑∂x²
        𝒑₀ᵀ𝗠⁻¹∂²𝒑∂y² = 𝒑₀ᵀ𝗠⁻¹*∂²𝒑∂y²
        𝒑₀ᵀ𝗠⁻¹∂²𝒑∂z² = 𝒑₀ᵀ𝗠⁻¹*∂²𝒑∂z²
        𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂x = 𝒑₀ᵀ∂𝗠⁻¹∂x*∂𝒑∂x
        𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂y = 𝒑₀ᵀ∂𝗠⁻¹∂y*∂𝒑∂y
        𝒑₀ᵀ∂𝗠⁻¹∂z∂𝒑∂z = 𝒑₀ᵀ∂𝗠⁻¹∂z*∂𝒑∂z
        𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂y = 𝒑₀ᵀ∂𝗠⁻¹∂x*∂𝒑∂y
        𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂x = 𝒑₀ᵀ∂𝗠⁻¹∂y*∂𝒑∂x
        𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂z = 𝒑₀ᵀ∂𝗠⁻¹∂x*∂𝒑∂z
        𝒑₀ᵀ∂𝗠⁻¹∂z∂𝒑∂x = 𝒑₀ᵀ∂𝗠⁻¹∂z*∂𝒑∂x
        𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂z = 𝒑₀ᵀ∂𝗠⁻¹∂y*∂𝒑∂z
        𝒑₀ᵀ∂𝗠⁻¹∂z∂𝒑∂y = 𝒑₀ᵀ∂𝗠⁻¹∂z*∂𝒑∂y
        𝒑₀ᵀ∂²𝗠⁻¹∂x∂y𝒑 = 𝒑₀ᵀ∂²𝗠⁻¹∂x∂y*𝒑
        𝒑₀ᵀ∂²𝗠⁻¹∂x∂z𝒑 = 𝒑₀ᵀ∂²𝗠⁻¹∂x∂z*𝒑
        𝒑₀ᵀ∂²𝗠⁻¹∂y∂z𝒑 = 𝒑₀ᵀ∂²𝗠⁻¹∂y∂z*𝒑
        𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x∂y = 𝒑₀ᵀ𝗠⁻¹*∂²𝒑∂x∂y
        𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x∂z = 𝒑₀ᵀ𝗠⁻¹*∂²𝒑∂x∂z
        𝒑₀ᵀ𝗠⁻¹∂²𝒑∂y∂z = 𝒑₀ᵀ𝗠⁻¹*∂²𝒑∂y∂z

        𝝭[i] = 𝒑₀ᵀ𝗠⁻¹𝒑*𝜙
        ∂𝝭∂x[i] = 𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂𝜙∂x
        ∂𝝭∂y[i] = 𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂𝜙∂y
        ∂𝝭∂z[i] = 𝒑₀ᵀ∂𝗠⁻¹∂z𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂z*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂𝜙∂z

        ∂²𝝭∂x²[i] = 𝒑₀ᵀ∂²𝗠⁻¹∂x²𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x²*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂²𝜙∂x² + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂x*𝜙 + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*∂𝜙∂x + 2.0*𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*∂𝜙∂x

        ∂²𝝭∂y²[i] = 𝒑₀ᵀ∂²𝗠⁻¹∂y²𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂y²*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂²𝜙∂y² + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂y*𝜙 + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*∂𝜙∂y + 2.0*𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*∂𝜙∂y

        ∂²𝝭∂z²[i] = 𝒑₀ᵀ∂²𝗠⁻¹∂z²𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂z²*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂²𝜙∂z² + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂z∂𝒑∂z*𝜙 + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂z𝒑*∂𝜙∂z + 2.0*𝒑₀ᵀ𝗠⁻¹∂𝒑∂z*∂𝜙∂z

        ∂²𝝭∂x∂y[i] = 𝒑₀ᵀ∂²𝗠⁻¹∂x∂y𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x∂y*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂²𝜙∂x∂y + 𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂y*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂x*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*∂𝜙∂y + 𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*∂𝜙∂x + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*∂𝜙∂y +𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*∂𝜙∂x

        ∂²𝝭∂x∂z[i] = 𝒑₀ᵀ∂²𝗠⁻¹∂x∂z*𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x∂z*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂²𝜙∂x∂z + 𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂z*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂z∂𝒑∂x*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*∂𝜙∂z + 𝒑₀ᵀ∂𝗠⁻¹∂z𝒑*∂𝜙∂x + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*∂𝜙∂z +𝒑₀ᵀ𝗠⁻¹∂𝒑∂z*∂𝜙∂x

        ∂²𝝭∂y∂z[i] = 𝒑₀ᵀ∂²𝗠⁻¹∂y∂z𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂y∂z*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂²𝜙∂y∂z + 𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂z*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂z∂𝒑∂y*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*∂𝜙∂z + 𝒑₀ᵀ∂𝗠⁻¹∂z𝒑*∂𝜙∂y + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*∂𝜙∂z +𝒑₀ᵀ𝗠⁻¹∂𝒑∂z*∂𝜙∂y
    end
    return 𝝭, ∂𝝭∂x, ∂𝝭∂y, ∂²𝝭∂x², ∂²𝝭∂x∂y, ∂²𝝭∂y², ∂𝝭∂z, ∂²𝝭∂z², ∂²𝝭∂x∂z, ∂²𝝭∂y∂z
end

function get∇³𝝭(ap::ReproducingKernel,𝒙::NTuple{3,Float64})
    𝓒 = ap.𝓒
    𝝭 = ap.𝝭[:∂1]
    ∂𝝭∂x = ap.𝝭[:∂x]
    ∂𝝭∂y = ap.𝝭[:∂y]
    ∂²𝝭∂x² = ap.𝝭[:∂x²]
    ∂²𝝭∂x∂y = ap.𝝭[:∂x∂y]
    ∂²𝝭∂y² = ap.𝝭[:∂y²]
    ∂³𝝭∂x³ = ap.𝝭[:∂x³]
    ∂³𝝭∂x²∂y = ap.𝝭[:∂x²∂y]
    ∂³𝝭∂x∂y² = ap.𝝭[:∂x∂y²]
    ∂³𝝭∂y³ = ap.𝝭[:∂y³]
    𝒑₀ᵀ𝗠⁻¹, 𝒑₀ᵀ∂𝗠⁻¹∂x, 𝒑₀ᵀ∂𝗠⁻¹∂y, 𝒑₀ᵀ∂²𝗠⁻¹∂x², 𝒑₀ᵀ∂²𝗠⁻¹∂x∂y, 𝒑₀ᵀ∂²𝗠⁻¹∂y², 𝒑₀ᵀ∂³𝗠⁻¹∂x³, 𝒑₀ᵀ∂³𝗠⁻¹∂x²∂y, 𝒑₀ᵀ∂³𝗠⁻¹∂x∂y², 𝒑₀ᵀ∂³𝗠⁻¹∂y³ = cal∇³𝗠!(ap,𝒙)
    for (i,𝒙ᵢ) in enumerate(𝓒)
        Δ𝒙 = 𝒙 - 𝒙ᵢ
        𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂²𝒑∂x², ∂²𝒑∂x∂y, ∂²𝒑∂y², ∂³𝒑∂x³, ∂³𝒑∂x²∂y, ∂³𝒑∂x∂y², ∂³𝒑∂y³ = get∇³𝒑(ap,Δ𝒙)
        𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂²𝜙∂x², ∂²𝜙∂x∂y, ∂²𝜙∂y², ∂³𝜙∂x³, ∂³𝜙∂x²∂y, ∂³𝜙∂x∂y², ∂³𝜙∂y³ = get∇³𝜙(ap,𝒙ᵢ,Δ𝒙)
        𝒑₀ᵀ𝗠⁻¹𝒑 = 𝒑₀ᵀ𝗠⁻¹*𝒑
        𝒑₀ᵀ∂𝗠⁻¹∂x𝒑 = 𝒑₀ᵀ∂𝗠⁻¹∂x*𝒑
        𝒑₀ᵀ∂𝗠⁻¹∂y𝒑 = 𝒑₀ᵀ∂𝗠⁻¹∂y*𝒑
        𝒑₀ᵀ𝗠⁻¹∂𝒑∂x = 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂x
        𝒑₀ᵀ𝗠⁻¹∂𝒑∂y = 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂y
        𝒑₀ᵀ∂²𝗠⁻¹∂x²𝒑 = 𝒑₀ᵀ∂²𝗠⁻¹∂x²*𝒑
        𝒑₀ᵀ∂²𝗠⁻¹∂x∂y𝒑 = 𝒑₀ᵀ∂²𝗠⁻¹∂x∂y*𝒑
        𝒑₀ᵀ∂²𝗠⁻¹∂y²𝒑 = 𝒑₀ᵀ∂²𝗠⁻¹∂y²*𝒑
        𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x² = 𝒑₀ᵀ𝗠⁻¹*∂²𝒑∂x²
        𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x∂y = 𝒑₀ᵀ𝗠⁻¹*∂²𝒑∂x∂y
        𝒑₀ᵀ𝗠⁻¹∂²𝒑∂y² = 𝒑₀ᵀ𝗠⁻¹*∂²𝒑∂y²
        𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂x = 𝒑₀ᵀ∂𝗠⁻¹∂x*∂𝒑∂x
        𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂y = 𝒑₀ᵀ∂𝗠⁻¹∂x*∂𝒑∂y
        𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂x = 𝒑₀ᵀ∂𝗠⁻¹∂y*∂𝒑∂x
        𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂y = 𝒑₀ᵀ∂𝗠⁻¹∂y*∂𝒑∂y
        𝒑₀ᵀ∂³𝗠⁻¹∂x³𝒑 = 𝒑₀ᵀ∂³𝗠⁻¹∂x³*𝒑
        𝒑₀ᵀ∂³𝗠⁻¹∂x²∂y𝒑 = 𝒑₀ᵀ∂³𝗠⁻¹∂x²∂y*𝒑
        𝒑₀ᵀ∂³𝗠⁻¹∂x∂y²𝒑 = 𝒑₀ᵀ∂³𝗠⁻¹∂x∂y²*𝒑
        𝒑₀ᵀ∂³𝗠⁻¹∂y³𝒑 = 𝒑₀ᵀ∂³𝗠⁻¹∂y³*𝒑
        𝒑₀ᵀ𝗠⁻¹∂³𝒑∂x³ = 𝒑₀ᵀ𝗠⁻¹*∂³𝒑∂x³
        𝒑₀ᵀ𝗠⁻¹∂³𝒑∂x²∂y = 𝒑₀ᵀ𝗠⁻¹*∂³𝒑∂x²∂y
        𝒑₀ᵀ𝗠⁻¹∂³𝒑∂x∂y² = 𝒑₀ᵀ𝗠⁻¹*∂³𝒑∂x∂y²
        𝒑₀ᵀ𝗠⁻¹∂³𝒑∂y³ = 𝒑₀ᵀ𝗠⁻¹*∂³𝒑∂y³
        𝒑₀ᵀ∂²𝗠⁻¹∂x²∂𝒑∂x = 𝒑₀ᵀ∂²𝗠⁻¹∂x²*∂𝒑∂x
        𝒑₀ᵀ∂𝗠⁻¹∂x∂²𝒑∂x² = 𝒑₀ᵀ∂𝗠⁻¹∂x*∂²𝒑∂x²
        𝒑₀ᵀ∂²𝗠⁻¹∂x²∂𝒑∂y = 𝒑₀ᵀ∂²𝗠⁻¹∂x²*∂𝒑∂y
        𝒑₀ᵀ∂𝗠⁻¹∂y∂²𝒑∂x² = 𝒑₀ᵀ∂𝗠⁻¹∂y*∂²𝒑∂x²
        𝒑₀ᵀ∂²𝗠⁻¹∂x∂y∂𝒑∂x = 𝒑₀ᵀ∂²𝗠⁻¹∂x∂y*∂𝒑∂x
        𝒑₀ᵀ∂𝗠⁻¹∂x∂²𝒑∂x∂y = 𝒑₀ᵀ∂𝗠⁻¹∂x*∂²𝒑∂x∂y
        𝒑₀ᵀ∂²𝗠⁻¹∂y²∂𝒑∂x = 𝒑₀ᵀ∂²𝗠⁻¹∂y²*∂𝒑∂x
        𝒑₀ᵀ∂𝗠⁻¹∂x∂²𝒑∂y² = 𝒑₀ᵀ∂𝗠⁻¹∂x*∂²𝒑∂y²
        𝒑₀ᵀ∂²𝗠⁻¹∂x∂y∂𝒑∂y = 𝒑₀ᵀ∂²𝗠⁻¹∂x∂y*∂𝒑∂y
        𝒑₀ᵀ∂𝗠⁻¹∂y∂²𝒑∂x∂y = 𝒑₀ᵀ∂𝗠⁻¹∂y*∂²𝒑∂x∂y
        𝒑₀ᵀ∂²𝗠⁻¹∂y²∂𝒑∂y = 𝒑₀ᵀ∂²𝗠⁻¹∂y²*∂𝒑∂y
        𝒑₀ᵀ∂𝗠⁻¹∂y∂²𝒑∂y² = 𝒑₀ᵀ∂𝗠⁻¹∂y*∂²𝒑∂y²

        𝝭[i] = 𝒑₀ᵀ𝗠⁻¹𝒑*𝜙
        ∂𝝭∂x[i] = 𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂𝜙∂x
        ∂𝝭∂y[i] = 𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂𝜙∂y

        ∂²𝝭∂x²[i] = 𝒑₀ᵀ∂²𝗠⁻¹∂x²𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x²*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂²𝜙∂x² + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂x*𝜙 + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*∂𝜙∂x + 2.0*𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*∂𝜙∂x

        ∂²𝝭∂y²[i] = 𝒑₀ᵀ∂²𝗠⁻¹∂y²𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂y²*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂²𝜙∂y² + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂y*𝜙 + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*∂𝜙∂y + 2.0*𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*∂𝜙∂y

        ∂²𝝭∂x∂y[i] = 𝒑₀ᵀ∂²𝗠⁻¹∂x∂y𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x∂y*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂²𝜙∂x∂y + 𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂y*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂x*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*∂𝜙∂y + 𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*∂𝜙∂x + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*∂𝜙∂y +𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*∂𝜙∂x

        ∂³𝝭∂x³[i] = 𝒑₀ᵀ∂³𝗠⁻¹∂x³𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂³𝒑∂x³*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂³𝜙∂x³ + 3*𝒑₀ᵀ∂²𝗠⁻¹∂x²∂𝒑∂x*𝜙  + 3*𝒑₀ᵀ∂𝗠⁻¹∂x∂²𝒑∂x²*𝜙 + 3*𝒑₀ᵀ∂²𝗠⁻¹∂x²𝒑*∂𝜙∂x + 3*𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*∂²𝜙∂x² + 3*𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x²*∂𝜙∂x + 3*𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*∂²𝜙∂x² + 6*𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂x*∂𝜙∂x

        ∂³𝝭∂x²∂y[i] = 𝒑₀ᵀ∂³𝗠⁻¹∂x²∂y𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂³𝒑∂x²∂y*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂³𝜙∂x²∂y + 2*𝒑₀ᵀ∂²𝗠⁻¹∂x∂y∂𝒑∂x*𝜙 + 𝒑₀ᵀ∂²𝗠⁻¹∂x²∂𝒑∂y*𝜙 + 2*𝒑₀ᵀ∂𝗠⁻¹∂x∂²𝒑∂x∂y*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂y∂²𝒑∂x²*𝜙 + 2*𝒑₀ᵀ∂²𝗠⁻¹∂x∂y𝒑*∂𝜙∂x + 𝒑₀ᵀ∂²𝗠⁻¹∂x²𝒑*∂𝜙∂y + 2*𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*∂²𝜙∂x∂y + 𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*∂²𝜙∂x² + 2*𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x∂y*∂𝜙∂x + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x²*∂𝜙∂y + 2*𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*∂²𝜙∂x∂y + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*∂²𝜙∂x² + 2*𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂x*∂𝜙∂x + 2*𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂y*∂𝜙∂x + 2*𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂x*∂𝜙∂y

        ∂³𝝭∂x∂y²[i] = 𝒑₀ᵀ∂³𝗠⁻¹∂x∂y²𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂³𝒑∂x∂y²*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂³𝜙∂x∂y² + 2*𝒑₀ᵀ∂²𝗠⁻¹∂x∂y∂𝒑∂y*𝜙 + 𝒑₀ᵀ∂²𝗠⁻¹∂y²∂𝒑∂x*𝜙 + 2*𝒑₀ᵀ∂𝗠⁻¹∂y∂²𝒑∂x∂y*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂x∂²𝒑∂y²*𝜙 + 2*𝒑₀ᵀ∂²𝗠⁻¹∂x∂y𝒑*∂𝜙∂y + 𝒑₀ᵀ∂²𝗠⁻¹∂y²𝒑*∂𝜙∂x + 2*𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*∂²𝜙∂x∂y + 𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*∂²𝜙∂y² + 2*𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x∂y*∂𝜙∂y + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂y²*∂𝜙∂x + 2*𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*∂²𝜙∂x∂y + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*∂²𝜙∂y² + 2*𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂y*∂𝜙∂y + 2*𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂x*∂𝜙∂y + 2*𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂y*∂𝜙∂x

        ∂³𝝭∂y³[i] = 𝒑₀ᵀ∂³𝗠⁻¹∂y³𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂³𝒑∂y³*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂³𝜙∂y³ + 3*𝒑₀ᵀ∂²𝗠⁻¹∂y²∂𝒑∂y*𝜙  + 3*𝒑₀ᵀ∂𝗠⁻¹∂y∂²𝒑∂y²*𝜙 + 3*𝒑₀ᵀ∂²𝗠⁻¹∂y²𝒑*∂𝜙∂y + 3*𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*∂²𝜙∂y² + 3*𝒑₀ᵀ𝗠⁻¹∂²𝒑∂y²*∂𝜙∂y + 3*𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*∂²𝜙∂y² + 6*𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂y*∂𝜙∂y
    end
    return 𝝭, ∂𝝭∂x, ∂𝝭∂y, ∂²𝝭∂x², ∂²𝝭∂x∂y, ∂²𝝭∂y², ∂³𝝭∂x³, ∂³𝝭∂x²∂y, ∂³𝝭∂x∂y², ∂³𝝭∂y³
end

function get𝝭(ap::ReproducingKernel{𝝃,𝒑̄,𝑠,𝜙̄,:Node},𝒙::NTuple{3,Float64},index::Vector{Int}) where {𝝃<:AbstractNode,𝒑̄,𝑠,𝜙̄}
    𝝭 = ap.𝝭[:∂1]
    𝒑₀ᵀ𝗠⁻¹= cal𝗠!(ap,𝒙)
    for i in 1:length(index)
        𝒙ᵢ = ap.𝓒[index[i]]
        Δ𝒙 = 𝒙 - 𝒙ᵢ
        𝒑 = get𝒑(ap,Δ𝒙)
        𝜙 = get𝜙(ap,𝒙ᵢ,Δ𝒙)
        𝝭[i] = 𝒑₀ᵀ𝗠⁻¹*𝒑*𝜙
    end
    return 𝝭
end

function get∇𝝭(ap::ReproducingKernel{𝝃,𝒑̄,𝑠,𝜙̄,:Node},𝒙::NTuple{3,Float64},index::Vector{Int}) where {𝝃<:AbstractNode,𝒑̄,𝑠,𝜙̄}
    𝝭 = ap.𝝭[:∂1]
    ∂𝝭∂x = ap.𝝭[:∂x]
    ∂𝝭∂y = ap.𝝭[:∂y]
    ∂𝝭∂z = ap.𝝭[:∂z]
    𝒑₀ᵀ𝗠⁻¹, 𝒑₀ᵀ∂𝗠⁻¹∂x, 𝒑₀ᵀ∂𝗠⁻¹∂y, 𝒑₀ᵀ∂𝗠⁻¹∂z= cal∇𝗠!(ap,𝒙)
    for i in 1:length(index)
        𝒙ᵢ = ap.𝓒[index[i]]
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
function set∇²𝝭!(aps::Vector{T}) where T<:ReproducingKernel{SNode}
    for ap in aps
        set∇²𝝭!(ap)
    end
end
function set∇³𝝭!(aps::Vector{T}) where T<:ReproducingKernel{SNode}
    for ap in aps
        set∇³𝝭!(ap)
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

function set∇²𝝭!(ap::ReproducingKernel{SNode})
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        i = ξ.id
        I = ξ.index[i]
        ξ̂ = Node(ξ)
        𝝭,∂𝝭∂x,∂𝝭∂y,∂²𝝭∂x²,∂²𝝭∂x∂y,∂²𝝭∂y²,∂𝝭∂z,∂²𝝭∂x∂z,∂²𝝭∂y∂z,∂²𝝭∂z² = get∇²𝝭(ap,ξ̂)
        for j in 1:length(𝓒)
            ξ.𝝭[:∂1][I+j] = 𝝭[j]
            ξ.𝝭[:∂x][I+j] = ∂𝝭∂x[j]
            ξ.𝝭[:∂y][I+j] = ∂𝝭∂y[j]
            ξ.𝝭[:∂x²][I+j] = ∂²𝝭∂x²[j]
            ξ.𝝭[:∂x∂y][I+j] = ∂²𝝭∂x∂y[j]
            ξ.𝝭[:∂y²][I+j] = ∂²𝝭∂y²[j]
            ξ.𝝭[:∂z][I+j] = ∂𝝭∂z[j]
            ξ.𝝭[:∂x∂z][I+j] = ∂²𝝭∂x∂z[j]
            ξ.𝝭[:∂y∂z][I+j] = ∂²𝝭∂y∂z[j]
        end
    end
end

function set∇³𝝭!(ap::ReproducingKernel{SNode})
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        i = ξ.id
        I = ξ.index[i]
        ξ̂ = Node(ξ)
        𝝭,∂𝝭∂x,∂𝝭∂y,∂²𝝭∂x²,∂²𝝭∂x∂y,∂²𝝭∂y²,∂³𝝭∂x³,∂³𝝭∂x²∂y,∂³𝝭∂x∂y²,∂³𝝭∂y³ = get∇³𝝭(ap,ξ̂)
        for j in 1:length(𝓒)
            ξ.𝝭[:∂1][I+j] = 𝝭[j]
            ξ.𝝭[:∂x][I+j] = ∂𝝭∂x[j]
            ξ.𝝭[:∂y][I+j] = ∂𝝭∂y[j]
            ξ.𝝭[:∂x²][I+j] = ∂²𝝭∂x²[j]
            ξ.𝝭[:∂x∂y][I+j] = ∂²𝝭∂x∂y[j]
            ξ.𝝭[:∂y²][I+j] = ∂²𝝭∂y²[j]
            ξ.𝝭[:∂x³][I+j] = ∂³𝝭∂x³[j]
            ξ.𝝭[:∂x²∂y][I+j] = ∂³𝝭∂x²∂y[j]
            ξ.𝝭[:∂x∂y²][I+j] = ∂³𝝭∂x∂y²[j]
            ξ.𝝭[:∂y³][I+j] = ∂³𝝭∂y³[j]
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
    ∂𝝭∂x = ap.𝝭[:∂x]
    ∂𝝭∂y = ap.𝝭[:∂y]
    ∂𝝭∂z = ap.𝝭[:∂z]
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

function get∇²𝝭(ap::ReproducingKernel,ξ::SNode)
    𝝭 = ap.𝝭[:∂1]
    ∂𝝭∂x = ap.𝝭[:∂x]
    ∂𝝭∂y = ap.𝝭[:∂y]
    ∂𝝭∂z = ap.𝝭[:∂z]
    ∂²𝝭∂x² = ap.𝝭[:∂x²]
    ∂²𝝭∂x∂y = ap.𝝭[:∂x∂y]
    ∂²𝝭∂y² = ap.𝝭[:∂y²]
    ∂²𝝭∂x∂z = ap.𝝭[:∂x∂z]
    ∂²𝝭∂y∂z = ap.𝝭[:∂y∂z]
    ∂²𝝭∂z² = ap.𝝭[:∂z²]
    i = ξ.id
    index = ξ.index
    for j in 1:length(ap.𝓒)
        𝝭[j] = ξ.𝝭[:∂1][index[i]+j]
        ∂𝝭∂x[j] = ξ.𝝭[:∂x][index[i]+j]
        ∂𝝭∂y[j] = ξ.𝝭[:∂y][index[i]+j]
        ∂𝝭∂z[j] = ξ.𝝭[:∂z][index[i]+j]
        ∂²𝝭∂x²[j] = ξ.𝝭[:∂x²][index[i]+j]
        ∂²𝝭∂x∂y[j] = ξ.𝝭[:∂x∂y][index[i]+j]
        ∂²𝝭∂y²[j] = ξ.𝝭[:∂y²][index[i]+j]
        ∂²𝝭∂x∂z[j] = ξ.𝝭[:∂x∂z][index[i]+j]
        ∂²𝝭∂y∂z[j] = ξ.𝝭[:∂y∂z][index[i]+j]
        ∂²𝝭∂z²[j] = ξ.𝝭[:∂z²][index[i]+j]
    end
    return 𝝭, ∂𝝭∂x, ∂𝝭∂y, ∂²𝝭∂x², ∂²𝝭∂x∂y, ∂²𝝭∂y², ∂𝝭∂z, ∂²𝝭∂x∂z, ∂²𝝭∂y∂z, ∂²𝝭∂z²
end

function get∇³𝝭(ap::ReproducingKernel,ξ::SNode)
    𝝭 = ap.𝝭[:∂1]
    ∂𝝭∂x = ap.𝝭[:∂x]
    ∂𝝭∂y = ap.𝝭[:∂y]
    ∂²𝝭∂x² = ap.𝝭[:∂x²]
    ∂²𝝭∂x∂y = ap.𝝭[:∂x∂y]
    ∂²𝝭∂y² = ap.𝝭[:∂y²]
    ∂³𝝭∂x³ = ap.𝝭[:∂x³]
    ∂³𝝭∂x²∂y = ap.𝝭[:∂x²∂y]
    ∂³𝝭∂x∂y² = ap.𝝭[:∂x∂y²]
    ∂³𝝭∂y³ = ap.𝝭[:∂y³]
    i = ξ.id
    index = ξ.index
    for j in 1:length(ap.𝓒)
        𝝭[j] = ξ.𝝭[:∂1][index[i]+j]
        ∂𝝭∂x[j] = ξ.𝝭[:∂x][index[i]+j]
        ∂𝝭∂y[j] = ξ.𝝭[:∂y][index[i]+j]
        ∂²𝝭∂x²[j] = ξ.𝝭[:∂x²][index[i]+j]
        ∂²𝝭∂x∂y[j] = ξ.𝝭[:∂x∂y][index[i]+j]
        ∂²𝝭∂y²[j] = ξ.𝝭[:∂y²][index[i]+j]
        ∂³𝝭∂x³[j] = ξ.𝝭[:∂x³][index[i]+j]
        ∂³𝝭∂x²∂y[j] = ξ.𝝭[:∂x²∂y][index[i]+j]
        ∂³𝝭∂x∂y²[j] = ξ.𝝭[:∂x∂y²][index[i]+j]
        ∂³𝝭∂y³[j] = ξ.𝝭[:∂y³][index[i]+j]
    end
    return 𝝭, ∂𝝭∂x, ∂𝝭∂y, ∂²𝝭∂x², ∂²𝝭∂x∂y, ∂²𝝭∂y², ∂³𝝭∂x³, ∂³𝝭∂x²∂y, ∂³𝝭∂x∂y², ∂³𝝭∂y³
end
function get∇̄𝝭(ap::ReproducingKernel,ξ::SNode)
    𝝭 = ap.𝝭[:∂1]
    ∂𝝭∂x = ap.𝝭[:∂x]
    ∂𝝭∂y = ap.𝝭[:∂y]
    ∂𝝭∂z = ap.𝝭[:∂z]
    i = ξ.id
    index = ξ.index
    for j in 1:length(ap.𝓒)
        𝝭[j] = ξ.𝝭[:∂1][index[i]+j]
        ∂𝝭∂x[j] = ξ.𝝭[:∂̄x][index[i]+j]
        ∂𝝭∂y[j] = haskey(ξ.𝝭,:∂̄y) ? ξ.𝝭[:∂̄y][index[i]+j] : 0.0
        ∂𝝭∂z[j] = haskey(ξ.𝝭,:∂̄z) ? ξ.𝝭[:∂̄z][index[i]+j] : 0.0
    end
    return 𝝭, ∂𝝭∂x, ∂𝝭∂y, ∂𝝭∂z
end
## RK gradient smoothing
function set∇̃𝝭!(aps::Vector{T}) where T<:ReproducingKernel
    for ap in aps
        set∇̃𝝭!(ap)
    end
end
set∇̃𝝭!(ap::T) where T<:ReproducingKernel{SNode} = set∇̃𝝭!(ap,ap)

function set∇̃²𝝭!(aps::Vector{T}) where T<:ReproducingKernel
    for ap in aps
        set∇̃²𝝭!(ap)
    end
end
set∇̃²𝝭!(ap::T) where T<:ReproducingKernel{SNode} = set∇̃²𝝭!(ap,ap)

function set∇̃𝝭!(gps::Vector{T},aps::Vector{S}) where {T<:ReproducingKernel,S<:ReproducingKernel}
    if length(gps) ≠ length(aps)
        error("Miss match element numbers")
    else
        for i in 1:length(gps)
            set∇̃𝝭!(gps[i],aps[i])
        end
    end
end

function set∇̃²𝝭!(gps::Vector{T},aps::Vector{S}) where {T<:ReproducingKernel,S<:ReproducingKernel}
    if length(gps) ≠ length(aps)
        error("Miss match element numbers")
    else
        for i in 1:length(gps)
            set∇̃²𝝭!(gps[i],aps[i])
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
function set∇̃𝝭!(gp::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Seg2},ap::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Seg2}) where {𝒑,𝑠,𝜙}
    n₁ =  1.0
    n₂ = -1.0
    𝓒 = gp.𝓒
    𝓖 = gp.𝓖
    for ξ̂ in 𝓖
        𝒒̂ = get𝒒(gp,ξ̂)
        𝗚⁻¹ = cal𝗚!(gp)
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
    n₁₂ = x₂-x₃;n₂₂ = x₃-x₁;n₃₂ = x₁-x₂
    𝓒 = gp.𝓒
    𝓖 = gp.𝓖
    for ξ̂ in 𝓖
        𝒒̂ = get𝒑₁(gp,ξ̂)
        𝗚⁻¹ = cal𝗚!(gp)
        𝒒̂ᵀ𝗚⁻¹ = 𝒒̂*𝗚⁻¹
        ∂𝝭∂x = gp.𝝭[:∂x]
        ∂𝝭∂y = gp.𝝭[:∂y]
        fill!(∂𝝭∂x,0.0)
        fill!(∂𝝭∂y,0.0)
        for ξ in ap.𝓖
            w = ξ.w
            wᵇ = ξ.wᵇ
            𝝭 = get𝝭(ap,ξ)
            𝒒, ∂𝒒∂ξ, ∂𝒒∂η = get∇𝒑₁(ap,ξ)
            𝒒̂ᵀ𝗚⁻¹𝒒 =  𝒒̂ᵀ𝗚⁻¹*𝒒
            𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ = 𝒒̂ᵀ𝗚⁻¹*∂𝒒∂ξ
            𝒒̂ᵀ𝗚⁻¹∂𝒒∂η = 𝒒̂ᵀ𝗚⁻¹*∂𝒒∂η
            nᵇ₁ = 0.0;nᵇ₂ = 0.0
            ξ.ξ == 0.0 ? (nᵇ₁ += n₁₁;nᵇ₂ += n₁₂) : nothing
            ξ.η == 0.0 ? (nᵇ₁ += n₂₁;nᵇ₂ += n₂₂) : nothing
            ξ.ξ+ξ.η ≈ 1.0 ? (nᵇ₁ += n₃₁;nᵇ₂ += n₃₂) : nothing
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

function set∇̃²𝝭!(gp::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Tri3},ap::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Tri3}) where {𝒑,𝑠,𝜙}
    x₁ = gp.𝓒[1].x;y₁ = gp.𝓒[1].y
    x₂ = gp.𝓒[2].x;y₂ = gp.𝓒[2].y
    x₃ = gp.𝓒[3].x;y₃ = gp.𝓒[3].y
    𝐴 = get𝐴(gp)
    n₁₁ = y₃-y₂;n₂₁ = y₁-y₃;n₃₁ = y₂-y₁
    n₁₂ = x₂-x₃;n₂₂ = x₃-x₁;n₃₂ = x₁-x₂
    s₁₁ = -n₁₂;s₂₁ = -n₂₂;s₃₁ = -n₃₂
    s₁₂ =  n₁₁;s₂₂ =  n₂₁;s₃₂ =  n₃₁
    𝓒 = gp.𝓒
    𝓖 = gp.𝓖
    for ξ̂ in 𝓖
        𝒒̂ = get𝒑₂(gp,ξ̂)
        𝗚⁻¹ = cal𝗚₂!(gp)
        𝒒̂ᵀ𝗚⁻¹ = 𝒒̂*𝗚⁻¹
        ∂²𝝭∂x² = gp.𝝭[:∂x²]
        ∂²𝝭∂x∂y = gp.𝝭[:∂x∂y]
        ∂²𝝭∂y² = gp.𝝭[:∂y²]
        fill!(∂²𝝭∂x²,0.0)
        fill!(∂²𝝭∂x∂y,0.0)
        fill!(∂²𝝭∂y²,0.0)
        for ξ in ap.𝓖
            w = ξ.w
            wᵇ = ξ.wᵇ
            𝝭,∂𝝭∂x,∂𝝭∂y = get∇𝝭(ap,ξ)
            𝒒, ∂𝒒∂ξ, ∂𝒒∂η, ∂²𝒒∂ξ², ∂²𝒒∂ξ∂η, ∂²𝒒∂η² = get∇²𝒑₂(ap,ξ)
            𝒒̂ᵀ𝗚⁻¹𝒒 =  𝒒̂ᵀ𝗚⁻¹*𝒒
            𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ = 𝒒̂ᵀ𝗚⁻¹*∂𝒒∂ξ
            𝒒̂ᵀ𝗚⁻¹∂𝒒∂η = 𝒒̂ᵀ𝗚⁻¹*∂𝒒∂η
            𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ² = 𝒒̂ᵀ𝗚⁻¹*∂²𝒒∂ξ²
            𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ∂η = 𝒒̂ᵀ𝗚⁻¹*∂²𝒒∂ξ∂η
            𝒒̂ᵀ𝗚⁻¹∂²𝒒∂η² = 𝒒̂ᵀ𝗚⁻¹*∂²𝒒∂η²

            q₁₁ = 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ²*n₁₁*n₁₁ + 2*𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ∂η*n₁₁*n₂₁ + 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂η²*n₂₁*n₂₁
            q₁₂ = 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ²*n₁₁*n₁₂ + 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ∂η*(n₁₁*n₂₂+n₁₂*n₂₁) + 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂η²*n₂₁*n₂₂
            q₂₂ = 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ²*n₁₂*n₁₂ + 2*𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ∂η*n₁₂*n₂₂ + 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂η²*n₂₂*n₂₂

            q₁ = 𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ*n₁₁ + 𝒒̂ᵀ𝗚⁻¹∂𝒒∂η*n₂₁
            q₂ = 𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ*n₁₂ + 𝒒̂ᵀ𝗚⁻¹∂𝒒∂η*n₂₂

            q₁n₁ = 0.0;q₂n₂ = 0.0;q₁n₂ = 0.0;q₂n₁ = 0.0
            mn₁₁n₁ = 0.0;mn₁₁n₂ = 0.0;mn₁₂n₁ = 0.0;mn₁₂n₂ = 0.0;mn₂₂n₁ = 0.0;mn₂₂n₂ = 0.0
            ms₁₁ = 0.0;ms₁₂ = 0.0;ms₂₂ = 0.0
            Δms₁₁ = 0.0;Δms₁₂ = 0.0;Δms₂₂ = 0.0
            𝐿₁² = n₁₁^2+n₁₂^2
            𝐿₂² = n₂₁^2+n₂₂^2
            𝐿₃² = n₃₁^2+n₃₂^2
            if ξ.ξ == 0.0
                q₁n₁ += q₁*n₁₁
                q₁n₂ += q₁*n₁₂
                q₂n₁ += q₂*n₁₁
                q₂n₂ += q₂*n₁₂
                mn₁₁n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₁₁*n₁₁*n₁₁/𝐿₁²
                mn₁₁n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₁₁*n₁₁*n₁₂/𝐿₁²
                mn₁₂n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₁₁*n₁₂*n₁₁/𝐿₁²
                mn₁₂n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₁₁*n₁₂*n₁₂/𝐿₁²
                mn₂₂n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₁₂*n₁₂*n₁₁/𝐿₁²
                mn₂₂n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₁₂*n₁₂*n₁₂/𝐿₁²
                ms₁₁ += (q₁*s₁₁+q₂*s₁₂)*n₁₁*s₁₁/𝐿₁²
                ms₁₂ += (q₁*s₁₁+q₂*s₁₂)*0.5*(n₁₁*s₁₂+n₁₂*s₁₁)/𝐿₁²
                ms₂₂ += (q₁*s₁₁+q₂*s₁₂)*n₁₂*s₁₂/𝐿₁²
                if ξ.η == 0.0
                    Δms₁₁ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₁₁*s₁₁/𝐿₁²-n₂₁*s₂₁/𝐿₂²)
                    Δms₂₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₁₂*s₁₂/𝐿₁²-n₂₂*s₂₂/𝐿₂²)
                    Δms₁₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*0.5*((n₁₁*s₁₂+n₁₂*s₁₁)/𝐿₁²-(n₂₁*s₂₂+n₂₂*s₂₁)/𝐿₂²)
                end
            end
            if  ξ.η == 0.0
                q₁n₁ += q₁*n₂₁
                q₁n₂ += q₁*n₂₂
                q₂n₁ += q₂*n₂₁
                q₂n₂ += q₂*n₂₂
                mn₁₁n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₂₁*n₂₁*n₂₁/𝐿₂²
                mn₁₁n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₂₁*n₂₁*n₂₂/𝐿₂²
                mn₁₂n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₂₁*n₂₂*n₂₁/𝐿₂²
                mn₁₂n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₂₁*n₂₂*n₂₂/𝐿₂²
                mn₂₂n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₂₂*n₂₂*n₂₁/𝐿₂²
                mn₂₂n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₂₂*n₂₂*n₂₂/𝐿₂²
                ms₁₁ += (q₁*s₂₁+q₂*s₂₂)*n₂₁*s₂₁/𝐿₂²
                ms₁₂ += (q₁*s₂₁+q₂*s₂₂)*0.5*(n₂₁*s₂₂+n₂₂*s₂₁)/𝐿₂²
                ms₂₂ += (q₁*s₂₁+q₂*s₂₂)*n₂₂*s₂₂/𝐿₂²
                if ξ.ξ+ξ.η ≈ 1.0
                    Δms₁₁ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₂₁*s₂₁/𝐿₂²-n₃₁*s₃₁/𝐿₃²)
                    Δms₂₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₂₂*s₂₂/𝐿₂²-n₃₂*s₃₂/𝐿₃²)
                    Δms₁₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*0.5*((n₂₁*s₂₂+n₂₂*s₂₁)/𝐿₂²-(n₃₁*s₃₂+n₃₂*s₃₁)/𝐿₃²)
                end
            end
            if ξ.ξ+ξ.η ≈ 1.0
                q₁n₁ += q₁*n₃₁
                q₁n₂ += q₁*n₃₂
                q₂n₁ += q₂*n₃₁
                q₂n₂ += q₂*n₃₂
                mn₁₁n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₃₁*n₃₁*n₃₁/𝐿₃²
                mn₁₁n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₃₁*n₃₁*n₃₂/𝐿₃²
                mn₁₂n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₃₁*n₃₂*n₃₁/𝐿₃²
                mn₁₂n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₃₁*n₃₂*n₃₂/𝐿₃²
                mn₂₂n₁ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₃₂*n₃₂*n₃₁/𝐿₃²
                mn₂₂n₂ += 𝒒̂ᵀ𝗚⁻¹𝒒*n₃₂*n₃₂*n₃₂/𝐿₃²
                ms₁₁ += (q₁*s₃₁+q₂*s₃₂)*n₃₁*s₃₁/𝐿₃²
                ms₁₂ += (q₁*s₃₁+q₂*s₃₂)*0.5*(n₃₁*s₃₂+n₃₂*s₃₁)/𝐿₃²
                ms₂₂ += (q₁*s₃₁+q₂*s₃₂)*n₃₂*s₃₂/𝐿₃²
                if ξ.ξ == 0.0
                    Δms₁₁ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₃₁*s₃₁/𝐿₃²-n₁₁*s₁₁/𝐿₁²)
                    Δms₂₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₃₂*s₃₂/𝐿₃²-n₁₂*s₁₂/𝐿₁²)
                    Δms₁₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*0.5*((n₃₁*s₃₂+n₃₂*s₃₁)/𝐿₃²-(n₁₁*s₁₂+n₁₂*s₁₁)/𝐿₁²)
                end
            end

            W₁₁₁ = mn₁₁n₁*wᵇ
            W₁₁₂ = mn₁₁n₂*wᵇ
            W₁₂₁ = mn₁₂n₁*wᵇ
            W₁₂₂ = mn₁₂n₂*wᵇ
            W₂₂₁ = mn₂₂n₁*wᵇ
            W₂₂₂ = mn₂₂n₂*wᵇ
            W₁₁ = (q₁₁*w + 2*(q₁n₁+ms₁₁)*wᵇ)/4/𝐴 + Δms₁₁
            W₁₂ = (q₁₂*w + (q₁n₂+q₂n₁+2*ms₁₂)*wᵇ)/4/𝐴 + Δms₁₂
            W₂₂ = (q₂₂*w + 2*(q₂n₂+ms₂₂)*wᵇ)/4/𝐴 + Δms₂₂
            for i in 1:length(𝓒)
                ∂²𝝭∂x²[i] += 𝝭[i]*W₁₁ + ∂𝝭∂x[i]*W₁₁₁ + ∂𝝭∂y[i]*W₁₁₂
                ∂²𝝭∂x∂y[i] += 𝝭[i]*W₁₂ + ∂𝝭∂x[i]*W₁₂₁ + ∂𝝭∂y[i]*W₁₂₂
                ∂²𝝭∂y²[i] += 𝝭[i]*W₂₂ + ∂𝝭∂x[i]*W₂₂₁ + ∂𝝭∂y[i]*W₂₂₂
            end
        end
        for i in 1:length(𝓒)
            ξ̂.𝝭[:∂x²][ξ̂.index[ξ̂.id]+i] = ∂²𝝭∂x²[i]
            ξ̂.𝝭[:∂x∂y][ξ̂.index[ξ̂.id]+i] = ∂²𝝭∂x∂y[i]
            ξ̂.𝝭[:∂y²][ξ̂.index[ξ̂.id]+i] = ∂²𝝭∂y²[i]
        end
    end
end
function set∇̃𝝭!(gp::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Tet4},ap::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Tet4}) where {𝒑,𝑠,𝜙}
    n₁₁
    𝗚⁻¹ = cal𝗚!(gp)
    𝓒 = gp.𝓒
    𝓖 = gp.𝓖
    for ξ̂ in 𝓖
        𝒒̂ = get𝒑₁(gp,ξ̂)
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
            𝒒, ∂𝒒∂ξ, ∂𝒒∂η, ∂𝒒∂γ = get∇𝒑₁(gp,ξ)
            b₁ = 𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ*n₁₁ + 𝒒̂ᵀ𝗚⁻¹∂𝒒∂η*n₂₁ + 𝒒̂ᵀ𝗚⁻¹∂𝒒∂γ*n₃₁
            b₂ = 𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ*n₁₂ + 𝒒̂ᵀ𝗚⁻¹∂𝒒∂η*n₂₂ + 𝒒̂ᵀ𝗚⁻¹∂𝒒∂γ*n₃₂
            b₂ = 𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ*n₁₃ + 𝒒̂ᵀ𝗚⁻¹∂𝒒∂η*n₂₃ + 𝒒̂ᵀ𝗚⁻¹∂𝒒∂γ*n₃₃
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

function set∇̄𝝭!(aps::Vector{T}) where T<:AbstractElement
    set𝝭!(aps)
    for ap in aps
        set∇̄𝝭!(ap)
    end
end

function set∇̄𝝭!(ap::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Seg2}) where {𝒑,𝑠,𝜙}
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ̂ in 𝓖
        𝒒̂ = get𝒑₁(ap,ξ̂)
        𝗚⁻¹ = cal𝗚!(ap)
        𝒒̂ᵀ𝗚⁻¹ = 𝒒̂*𝗚⁻¹
        ∂𝝭∂x = ap.𝝭[:∂x]
        fill!(∂𝝭∂x,0.0)
        for ξ in ap.𝓖
            w = ξ.w
            n = ξ.n₁
            𝝭 = get𝝭(ap,ξ)
            𝒒 = get𝒑₁(ap,ξ)
            W₁ = 𝒒̂ᵀ𝗚⁻¹*𝒒*n*w
            for i in 1:length(𝓒)
                ∂𝝭∂x[i] += 𝝭[i]*W₁
            end
        end
        for i in 1:length(𝓒)
            ξ̂.𝝭[:∂̄x][ξ̂.index[ξ̂.id]+i] = ∂𝝭∂x[i]
        end
    end
end

function set∇̄𝝭!(ap::ReproducingKernel{SNode,𝒑,𝑠,𝜙,:Tri3}) where {𝒑,𝑠,𝜙}
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ̂ in 𝓖
        𝒒̂ = get𝒑₁(ap,ξ̂)
        𝗚⁻¹ = cal𝗚!(ap)
        𝒒̂ᵀ𝗚⁻¹ = 𝒒̂*𝗚⁻¹
        ∂𝝭∂x = ap.𝝭[:∂x]
        ∂𝝭∂y = ap.𝝭[:∂y]
        fill!(∂𝝭∂x,0.0)
        fill!(∂𝝭∂y,0.0)
        for ξ in ap.𝓖
            w = ξ.w
            n₁ = ξ.n₁
            n₂ = ξ.n₂
            𝝭 = get𝝭(ap,ξ)
            𝒒 = get𝒑₁(ap,ξ)
            𝒒̂ᵀ𝗚⁻¹𝒒 = 𝒒̂ᵀ𝗚⁻¹*𝒒
            W₁ = 𝒒̂ᵀ𝗚⁻¹𝒒*n₁*w
            W₂ = 𝒒̂ᵀ𝗚⁻¹𝒒*n₂*w
            for i in 1:length(𝓒)
                ∂𝝭∂x[i] += 𝝭[i]*W₁
                ∂𝝭∂y[i] += 𝝭[i]*W₂
            end
        end
        for i in 1:length(𝓒)
            ξ̂.𝝭[:∂̄x][ξ̂.index[ξ̂.id]+i] = ∂𝝭∂x[i]
            ξ̂.𝝭[:∂̄y][ξ̂.index[ξ̂.id]+i] = ∂𝝭∂y[i]
        end
    end
end

@inline function set∇̃𝝭!(a::ReproducingKernel{SNode,𝒑,𝑠,𝜙,T},b::ReproducingKernel{SNode,𝒑,𝑠,𝜙,T},c::ReproducingKernel{SNode,𝒑,𝑠,𝜙,T}) where {𝒑,𝑠,𝜙,T}
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
function ReproducingKernel{𝝃,𝒑,𝑠,𝜙,T}(a::A,b::B,sp::Union{Nothing,SpatialPartition}=nothing) where {A<:ReproducingKernel,B<:ReproducingKernel,𝝃<:AbstractNode,𝒑,𝑠,𝜙,T}
    𝓒 = a.𝓒
    𝓖 = get𝓖(a,b)
    if 𝓖 ≠ nothing
        if 𝝃 == SNode
            n = length(a.𝓒)-length(b.𝓒)
            nₜ = length(𝓖)*n
            index = 𝓖[1].index
            𝝭 = 𝓖[1].𝝭
            for s in keys(𝝭)
                append!(𝝭[s],zeros(nₜ))
            end
            for ξ in 𝓖
                for i in 1:length(index)-ξ.id
                    index[ξ.id+i] += n
                end
            end
        end
        𝗠 = a.𝗠
        𝝭 = a.𝝭
        ap = ReproducingKernel{𝝃,𝒑,𝑠,𝜙,T}(𝓒,𝓖,𝗠,𝝭)
        sp ≠ nothing ? sp(ap) : nothing
        return ap
    else
        return nothing
    end
end

function ReproducingKernel{𝝃,𝒑,𝑠,𝜙,T}(as::Vector{A},bs::Vector{B}) where {𝝃<:AbstractNode,𝒑,𝑠,𝜙,T,A<:ReproducingKernel,B<:ReproducingKernel}
    aps = ReproducingKernel{𝝃,𝒑,𝑠,𝜙,T}[]
    for b in bs
        for a in as
            ap = ReproducingKernel{𝝃,𝒑,𝑠,𝜙,T}(a,b)
            ap ≠ nothing ? push!(aps,ap) : nothing
        end
    end
    return aps
end

## get∇𝑢
function get∇𝑢(ap::T,𝒙::NTuple{3,Float64},sp::S) where {T<:ReproducingKernel,S<:SpatialPartition}
    index = [sp(𝒙...)...]
    N,B₁,B₂,B₃ = get∇𝝭(ap,𝒙,index)
    u = 0.0
    ∂u∂x = 0.0
    ∂u∂y = 0.0
    ∂u∂z = 0.0
    for i in 1:length(index)
        id = index[i]
        x = ap.𝓒[id]
        u += N[i]*x.d
        ∂u∂x += B₁[i]*x.d
        ∂u∂y += B₂[i]*x.d
        ∂u∂z += B₃[i]*x.d
    end
    return u,∂u∂x,∂u∂y,∂u∂z
end

function get𝝐(ap::T,𝒙::NTuple{3,Float64},sp::S) where {T<:ReproducingKernel,S<:SpatialPartition}
    index = [sp(𝒙...)...]
    N,B₁,B₂ = get∇𝝭(ap,𝒙,index)
    u = 0.0
    ε₁₁ = 0.0
    ε₂₂ = 0.0
    ε₁₂ = 0.0
    for i in 1:length(index)
        id = index[i]
        x = ap.𝓒[id]
        u += N[i]*x.d
        ε₁₁ += B₁[i]*x.d₁
        ε₂₂ += B₂[i]*x.d₂
        ε₁₂ += B₁[i]*x.d₂ + B₂[i]*x.d₁
    end
    return u,ε₁₁,ε₂₂,ε₁₂
end
