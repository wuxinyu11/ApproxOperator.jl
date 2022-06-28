
struct SymMat
    n::Int
    m::Vector{Float64}
end
SymMat(n::Int) = SymMat(n,zeros(Int(n*(n+1)/2)))

@inline function getindex(A::SymMat,i::Int,j::Int)
    i > j ? A.m[Int(j+i*(i-1)/2)] : A.m[Int(i+j*(j-1)/2)]
end

@inline function setindex!(A::SymMat,val::Float64,i::Int,j::Int)
    i > j ? A.m[Int(j+i*(i-1)/2)] = val : A.m[Int(i+j*(j-1)/2)] = val
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
        for j in 1:i
            A[i,j] = sum(A[i,k]*A[k,j] for k in i:n)
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
    return A
end

function UAUᵀ!(A::SymMat,U::SymMat)
    n = A.n
    for i in 1:n
        for j in i:n
            A[i,j] = sum(U[i,k]*A[k,l]*U[j,l] for k in i:n for l in j:n)
        end
    end
    return A
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
    return A
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
    (sp::t)(x::T) where T = sp(x.x,x.y,x.z)
    function (sp::t)(xs::T...) where T
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
        T<:ReproducingKernel ? set_memory_𝝭!(aps) : nothing
    end
    function (sp::t)(apss::Any...)
        for aps in apss
            sp(aps)
        end
    end
end

```
ReproducingKernel
```
struct ReproducingKernel{𝑝,𝑠,𝜙,T,N₁,N₂}<:AbstractElement{T}
    𝓒::Vector{Node{N₁}}
    𝓖::Vector{Node{N₂}}
    𝗠::Dict{Symbol,SymMat}
end

## Basis Function
@inline get∇₁𝒑(ap::ReproducingKernel,x::Any) = get𝒑(ap,x), get∂𝒑∂x(ap,x)
@inline get∇₂𝒑(ap::ReproducingKernel,x::Any) = get𝒑(ap,x), get∂𝒑∂x(ap,x), get∂𝒑∂y(ap,x)
@inline get∇₃𝒑(ap::ReproducingKernel,x::Any) = get𝒑(ap,x), get∂𝒑∂x(ap,x), get∂𝒑∂y(ap,x), get∂𝒑∂z(ap,x)
@inline get∇²₁𝒑(ap::ReproducingKernel,x::Any) = get𝒑(ap,x), get∂𝒑∂x(ap,x), get∂²𝒑∂x²(ap,x)
@inline get∇²₂𝒑(ap::ReproducingKernel,x::Any) = get𝒑(ap,x), get∂𝒑∂x(ap,x), get∂𝒑∂y(ap,x), get∂²𝒑∂x²(ap,x), get∂²𝒑∂x∂y(ap,x), get∂²𝒑∂y²(ap,x)
@inline get∇²₃𝒑(ap::ReproducingKernel,x::Any) = get𝒑(ap,x), get∂𝒑∂x(ap,x), get∂𝒑∂y(ap,x), get∂²𝒑∂x²(ap,x), get∂²𝒑∂x∂y(ap,x), get∂²𝒑∂y²(ap,x), get∂𝒑∂z(ap,x), get∂²𝒑∂x∂z(ap,x), get∂²𝒑∂y∂z(ap,x), get∂²𝒑∂z²(ap,x)
@inline get∇³₁𝒑(ap::ReproducingKernel,x::Any) = get𝒑(ap,x), get∂𝒑∂x(ap,x), get∂²𝒑∂x²(ap,x), get∂³𝒑∂x³(ap,x)
@inline get∇³₂𝒑(ap::ReproducingKernel,x::Any) = get𝒑(ap,x), get∂𝒑∂x(ap,x), get∂𝒑∂y(ap,x), get∂²𝒑∂x²(ap,x), get∂²𝒑∂x∂y(ap,x), get∂²𝒑∂y²(ap,x), get∂³𝒑∂x³(ap,x), get∂³𝒑∂x²∂y(ap,x), get∂³𝒑∂x∂y²(ap,x), get∂³𝒑∂y³(ap,x)
@inline get∇∇²𝒑(ap::ReproducingKernel,x::Any) = get𝒑(ap,x), get∂³𝒑∂x³(ap,x), get∂³𝒑∂x²∂y(ap,x), get∂³𝒑∂x²∂y(ap,x), get∂³𝒑∂x∂y²(ap,x), get∂³𝒑∂x∂y²(ap,x), get∂³𝒑∂y³(ap,x)
@inline get∇𝒑₁(ap::ReproducingKernel{𝒑,𝑠,𝜙,:Seg2},ξ::Any) where {𝒑,𝑠,𝜙} = get𝒑₁(ap,ξ), get∂𝒑₁∂ξ(ap,ξ)
@inline get∇𝒑₁(ap::ReproducingKernel{𝒑,𝑠,𝜙,:Tri3},ξ::Any) where {𝒑,𝑠,𝜙} = get𝒑₁(ap,ξ), get∂𝒑₁∂ξ(ap,ξ), get∂𝒑₁∂η(ap,ξ)
@inline get∇𝒑₂(ap::ReproducingKernel{𝒑,𝑠,𝜙,:Tri3},ξ::Any) where {𝒑,𝑠,𝜙} = get𝒑₂(ap,ξ), get∂𝒑₂∂ξ(ap,ξ), get∂𝒑₂∂η(ap,ξ)
@inline get∇²𝒑₂(ap::ReproducingKernel{𝒑,𝑠,𝜙,:Tri3},ξ::Any) where {𝒑,𝑠,𝜙} = get𝒑₂(ap,ξ), get∂𝒑₂∂ξ(ap,ξ), get∂𝒑₂∂η(ap,ξ), get∂²𝒑₂∂ξ²(ap,ξ), get∂²𝒑₂∂ξ∂η(ap,ξ), get∂²𝒑₂∂η²(ap,ξ)

# ------------ Linear1D ---------------
@inline get𝑛𝒑(::ReproducingKernel{:Linear1D}) = 2
@inline get𝒑(::ReproducingKernel{:Linear1D},x::NTuple{3,Float64}) = (1.,x[1])
@inline get∂𝒑∂x(::ReproducingKernel{:Linear1D},::NTuple{3,Float64}) = (0.,1.)
@inline get∂𝒑∂y(::ReproducingKernel{:Linear1D},::Any) = (0.,0.)
@inline get∂𝒑∂z(::ReproducingKernel{:Linear1D},::Any) = (0.,0.)
@inline get∂²𝒑∂x²(::ReproducingKernel{:Linear1D},::Any) = (0.,0.)
@inline get∂²𝒑∂y²(::ReproducingKernel{:Linear1D},::Any) = (0.,0.)
@inline get∂²𝒑∂z²(::ReproducingKernel{:Linear1D},::Any) = (0.,0.)
@inline get∂²𝒑∂x∂y(::ReproducingKernel{:Linear1D},::Any) = (0.,0.)
@inline get∂²𝒑∂x∂z(::ReproducingKernel{:Linear1D},::Any) = (0.,0.)
@inline get∂²𝒑∂y∂z(::ReproducingKernel{:Linear1D},::Any) = (0.,0.)

@inline get𝑛𝒑₁(::ReproducingKernel{:Linear1D}) = 1
@inline get𝒑₁(::ReproducingKernel{:Linear1D},::Any) = (1.0,)
@inline get∂𝒑₁∂ξ(::ReproducingKernel{:Linear1D},::Any) = (0.0,)

# ------------ Quadaratic1D ---------------
@inline get𝑛𝒑(::ReproducingKernel{:Quadratic1D}) = 3
@inline get𝒑(::ReproducingKernel{:Quadratic1D},x::NTuple{3,Float64}) = (1.,x[1],x[1]^2)
@inline get∂𝒑∂x(::ReproducingKernel{:Quadratic1D},x::NTuple{3,Float64}) = (0.,1.,2*x[1])
@inline get∂𝒑∂y(::ReproducingKernel{:Quadratic1D},::Any) = (0.,0.,0.)
@inline get∂𝒑∂z(::ReproducingKernel{:Quadratic1D},::Any) = (0.,0.,0.)
@inline get∂²𝒑∂x²(::ReproducingKernel{:Quadratic1D},::Any) =(0.,0.,2.)
@inline get∂²𝒑∂y²(::ReproducingKernel{:Quadratic1D},::Any) =(0.,0.,0.)
@inline get∂²𝒑∂z²(::ReproducingKernel{:Quadratic1D},::Any) =(0.,0.,0.)
@inline get∂²𝒑∂x∂y(::ReproducingKernel{:Quadratic1D},::Any) =(0.,0.,0.)
@inline get∂²𝒑∂x∂z(::ReproducingKernel{:Quadratic1D},::Any) =(0.,0.,0.)
@inline get∂²𝒑∂y∂z(::ReproducingKernel{:Quadratic1D},::Any) =(0.,0.,0.)
@inline get∂³𝒑∂x³(::ReproducingKernel{:Quadratic1D},::Any) =(0.,0.,0.)
@inline get∂³𝒑∂x²∂y(::ReproducingKernel{:Quadratic1D},::Any) =(0.,0.,0.)
@inline get∂³𝒑∂x∂y²(::ReproducingKernel{:Quadratic1D},::Any) =(0.,0.,0.)
@inline get∂³𝒑∂y³(::ReproducingKernel{:Quadratic1D},::Any) =(0.,0.,0.)

@inline get𝑛𝒑₁(::ReproducingKernel{:Quadratic1D}) = 2
@inline get𝒑₁(ap::ReproducingKernel{:Quadratic1D},ξ::𝝃) = get𝒑₁(ap,ξ.ξ)
@inline get𝒑₁(::ReproducingKernel{:Quadratic1D},ξ::Float64) = (1.0,0.5*(1.0-ξ))
@inline get∂𝒑₁∂ξ(::ReproducingKernel{:Quadratic1D},::Any) = (0.0,1.0)

# ------------ Cubic1D ---------------
@inline get𝑛𝒑(::ReproducingKernel{:Cubic1D}) = 4
@inline get𝒑(::ReproducingKernel{:Cubic1D},x::NTuple{3,Float64}) = (1.,x[1],x[1]^2,x[1]^3)
@inline get∂𝒑∂x(::ReproducingKernel{:Cubic1D},x::NTuple{3,Float64}) = (0.,1.,2*x[1],3*x[1]^2)
@inline get∂𝒑∂y(::ReproducingKernel{:Cubic1D}, ::Any) = (0.,0.,0.,0.)
@inline get∂𝒑∂z(::ReproducingKernel{:Cubic1D}, ::Any) = (0.,0.,0.,0.)
@inline get∂²𝒑∂x²(::ReproducingKernel{:Cubic1D},x::NTuple{3,Float64}) = (0.,0.,2.,6*x[1])
@inline get∂²𝒑∂y²(::ReproducingKernel{:Cubic1D},x::NTuple{3,Float64}) = (0.,0.,0.,0.)
@inline get∂²𝒑∂z²(::ReproducingKernel{:Cubic1D},x::NTuple{3,Float64}) = (0.,0.,0.,0.)
@inline get∂²𝒑∂x∂y(::ReproducingKernel{:Cubic1D},x::NTuple{3,Float64}) = (0.,0.,0.,0.)
@inline get∂²𝒑∂x∂z(::ReproducingKernel{:Cubic1D},x::NTuple{3,Float64}) = (0.,0.,0.,0.)
@inline get∂²𝒑∂y∂z(::ReproducingKernel{:Cubic1D},x::NTuple{3,Float64}) = (0.,0.,0.,0.)
@inline get∂³𝒑∂x³(::ReproducingKernel{:Cubic1D},x::NTuple{3,Float64}) = (0.,0.,0.,6.)
@inline get∂³𝒑∂x²∂y(::ReproducingKernel{:Cubic1D},x::NTuple{3,Float64}) = (0.,0.,0.,0.)
@inline get∂³𝒑∂x∂y²(::ReproducingKernel{:Cubic1D},x::NTuple{3,Float64}) = (0.,0.,0.,0.)
@inline get∂³𝒑∂y³(::ReproducingKernel{:Cubic1D},x::NTuple{3,Float64}) = (0.,0.,0.,0.)

@inline get𝑛𝒑₁(::ReproducingKernel{:Cubic1D}) = 3
@inline get𝒑₁(ap::ReproducingKernel{:Cubic1D},ξ::𝝃) = get𝒑₁(ap,ξ.ξ)
@inline get𝒑₁(::ReproducingKernel{:Cubic1D},ξ::Float64) = (1.0,0.5*(1.0-ξ),0.25*(1.0-ξ)^2)
@inline get∂𝒑₁∂ξ(ap::ReproducingKernel{:Cubic1D},ξ::𝝃) = get∂𝒑₁∂ξ(ap,ξ.ξ)
@inline get∂𝒑₁∂ξ(::ReproducingKernel{:Cubic1D},ξ::Float64) = (0.,1.0,(1.0-ξ))

# ------------ Linear2D ---------------
@inline get𝑛𝒑(::ReproducingKernel{:Linear2D}) = 3
@inline get𝒑(::ReproducingKernel{:Linear2D},x::NTuple{3,Float64}) = (1.,x[1],x[2])
@inline get∂𝒑∂x(::ReproducingKernel{:Linear2D}, ::Any) = (0.,1.,0.)
@inline get∂𝒑∂y(::ReproducingKernel{:Linear2D}, ::Any) = (0.,0.,1.)
@inline get∂𝒑∂z(::ReproducingKernel{:Linear2D}, ::Any) = (0.,0.,0.)

@inline get𝑛𝒑₁(::ReproducingKernel{:Linear2D}) = 1
@inline get𝒑₁(ap::ReproducingKernel{:Linear2D},ξ::𝝃) = get𝒑₁(ap,ξ.ξ,ξ.η)
@inline get𝒑₁(::ReproducingKernel{:Linear2D},::Any,::Any) = (1.,)
@inline get∂𝒑₁∂ξ(ap::ReproducingKernel{:Linear2D},ξ::𝝃) = get∂𝒑₁∂ξ(ap,ξ.ξ,ξ.η)
@inline get∂𝒑₁∂ξ(::ReproducingKernel{:Linear2D},::Any,::Any) = (0.,)
@inline get∂𝒑₁∂η(ap::ReproducingKernel{:Linear2D},ξ::𝝃) = get∂𝒑₁∂η(ap,ξ.ξ,ξ.η)
@inline get∂𝒑₁∂η(::ReproducingKernel{:Linear2D},::Any,::Any) = (0.,)

# ------------ Quadratic2D ---------------
@inline get𝑛𝒑(::ReproducingKernel{:Quadratic2D}) = 6
@inline get𝒑(::ReproducingKernel{:Quadratic2D},x::NTuple{3,Float64}) = (1.,x[1],x[2],x[1]^2,x[1]*x[2],x[2]^2)
@inline get∂𝒑∂x(::ReproducingKernel{:Quadratic2D},x::NTuple{3,Float64}) = (0.,1.,0.,2*x[1],x[2],0.)
@inline get∂𝒑∂y(::ReproducingKernel{:Quadratic2D},x::NTuple{3,Float64}) = (0.,0.,1.,0.,x[1],2*x[2])
@inline get∂𝒑∂z(::ReproducingKernel{:Quadratic2D}, ::Any) = (0.,0.,0.,0.,0.,0.)
@inline get∂²𝒑∂x²(::ReproducingKernel{:Quadratic2D},x::NTuple{3,Float64}) = (0.,0.,0.,2.,0.,0.)
@inline get∂²𝒑∂y²(::ReproducingKernel{:Quadratic2D},x::NTuple{3,Float64}) = (0.,0.,0.,0.,0.,2.)
@inline get∂²𝒑∂z²(::ReproducingKernel{:Quadratic2D},x::NTuple{3,Float64}) = (0.,0.,0.,0.,0.,0.)
@inline get∂²𝒑∂x∂y(::ReproducingKernel{:Quadratic2D},x::NTuple{3,Float64}) = (0.,0.,0.,0.,1.,0.)
@inline get∂²𝒑∂x∂z(::ReproducingKernel{:Quadratic2D},x::NTuple{3,Float64}) = (0.,0.,0.,0.,0.,0.)
@inline get∂²𝒑∂y∂z(::ReproducingKernel{:Quadratic2D},x::NTuple{3,Float64}) = (0.,0.,0.,0.,0.,0.)
@inline get∂³𝒑∂x³(::ReproducingKernel{:Quadratic2D},x::NTuple{3,Float64}) = (0.,0.,0.,0.,0.,0.)
@inline get∂³𝒑∂x²∂y(::ReproducingKernel{:Quadratic2D},x::NTuple{3,Float64}) = (0.,0.,0.,0.,0.,0.)
@inline get∂³𝒑∂x∂y²(::ReproducingKernel{:Quadratic2D},x::NTuple{3,Float64}) = (0.,0.,0.,0.,0.,0.)
@inline get∂³𝒑∂y³(::ReproducingKernel{:Quadratic2D},x::NTuple{3,Float64}) = (0.,0.,0.,0.,0.,0.)

@inline get𝑛𝒑₁(::ReproducingKernel{:Quadratic2D}) = 3
@inline get𝒑₁(ap::ReproducingKernel{:Quadratic2D},ξ::𝝃) = get𝒑₁(ap,ξ.ξ,ξ.η)
@inline get𝒑₁(::ReproducingKernel{:Quadratic2D},ξ::Float64,η::Float64) = (1.,ξ,η)
@inline get∂𝒑₁∂ξ(::ReproducingKernel{:Quadratic2D},::Any) = (0.,1.,0.)
@inline get∂𝒑₁∂η(::ReproducingKernel{:Quadratic2D},::Any) = (0.,0.,1.)
@inline get∂²𝒑₁∂ξ²(::ReproducingKernel{:Quadratic2D},::Any) = (0.,0.,0.)
@inline get∂²𝒑₁∂ξ∂η(::ReproducingKernel{:Quadratic2D},::Any) = (0.,0.,0.)
@inline get∂²𝒑₁∂η²(::ReproducingKernel{:Quadratic2D},::Any) = (0.,0.,0.)
@inline get∂𝒑₁∂ξ(::ReproducingKernel{:Quadratic2D},::Any,::Any) = (0.,1.,0.)
@inline get∂𝒑₁∂η(::ReproducingKernel{:Quadratic2D},::Any,::Any) = (0.,0.,1.)

@inline get𝑛𝒑₂(::ReproducingKernel{:Quadratic2D}) = 1
@inline get𝒑₂(::ReproducingKernel{:Quadratic2D},::Any) = (1.,)
@inline get∂𝒑₂∂ξ(::ReproducingKernel{:Quadratic2D},::Any) = (0.,)
@inline get∂𝒑₂∂η(::ReproducingKernel{:Quadratic2D},::Any) = (0.,)
@inline get∂²𝒑₂∂ξ²(::ReproducingKernel{:Quadratic2D},::Any) = (0.,)
@inline get∂²𝒑₂∂ξ∂η(::ReproducingKernel{:Quadratic2D},::Any) = (0.,)

# ------------ Cubic2D ---------------
@inline get𝑛𝒑(::ReproducingKernel{:Cubic2D}) = 10
@inline get𝒑(::ReproducingKernel{:Cubic2D},x::NTuple{3,Float64}) =
(
    1., x[1], x[2], x[1]^2, x[1]*x[2], x[2]^2, x[1]^3, x[1]^2*x[2], x[1]*x[2]^2, x[2]^3
)
@inline get∂𝒑∂x(::ReproducingKernel{:Cubic2D},x::NTuple{3,Float64}) =
(
    0., 1., 0., 2*x[1], x[2], 0., 3*x[1]^2, 2*x[1]*x[2], x[2]^2, 0.
)
@inline get∂𝒑∂y(::ReproducingKernel{:Cubic2D},x::NTuple{3,Float64}) =
(
    0., 0., 1., 0., x[1], 2*x[2], 0., x[1]^2, 2*x[1]*x[2], 3*x[2]^2
)
@inline get∂²𝒑∂x²(::ReproducingKernel{:Cubic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 2., 0., 0., 6*x[1], 2*x[2], 0., 0.
)
@inline get∂²𝒑∂x∂y(::ReproducingKernel{:Cubic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 1., 0., 0., 2*x[1], 2*x[2], 0.
)
@inline get∂²𝒑∂y²(::ReproducingKernel{:Cubic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 2., 0., 0., 2*x[1], 6*x[2]
)
@inline get∂𝒑∂z(::ReproducingKernel{:Cubic2D},::Any) =
(
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
)
@inline get∂²𝒑∂x∂z(::ReproducingKernel{:Cubic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
)
@inline get∂²𝒑∂y∂z(::ReproducingKernel{:Cubic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
)
@inline get∂²𝒑∂z²(::ReproducingKernel{:Cubic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
)
@inline get∂³𝒑∂x³(::ReproducingKernel{:Cubic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 0., 6., 0., 0., 0.
)
@inline get∂³𝒑∂x²∂y(::ReproducingKernel{:Cubic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 0., 0., 2., 0., 0.
)
@inline get∂³𝒑∂x∂y²(::ReproducingKernel{:Cubic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 0., 0., 0., 2., 0.
)
@inline get∂³𝒑∂y³(::ReproducingKernel{:Cubic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 0., 0., 0., 0., 6.
)

@inline get𝑛𝒑₁(::ReproducingKernel{:Cubic2D}) = 6
@inline get𝒑₁(ap::ReproducingKernel{:Cubic2D},ξ::𝝃) = get𝒑₁(ap,ξ.ξ,ξ.η)
@inline get𝒑₁(::ReproducingKernel{:Cubic2D},ξ::Float64,η::Float64) = (1.,ξ,η,ξ^2,ξ*η,η^2)
@inline get∂𝒑₁∂ξ(ap::ReproducingKernel{:Cubic2D},ξ::𝝃) = get∂𝒑₁∂ξ(ap,ξ.ξ,ξ.η)
@inline get∂𝒑₁∂ξ(::ReproducingKernel{:Cubic2D},ξ::Float64,η::Float64) = (0.,1.,0.,2.0*ξ,η,0.)
@inline get∂𝒑₁∂η(ap::ReproducingKernel{:Cubic2D},ξ::𝝃) = get∂𝒑₁∂η(ap,ξ.ξ,ξ.η)
@inline get∂𝒑₁∂η(::ReproducingKernel{:Cubic2D},ξ::Float64,η::Float64) = (0.,0.,1.,0.,ξ,2.0*η)

@inline get𝑛𝒑₂(::ReproducingKernel{:Cubic2D}) = 3
@inline get𝒑₂(ap::ReproducingKernel{:Cubic2D},ξ::𝝃) = get𝒑₂(ap,ξ.ξ,ξ.η)
@inline get𝒑₂(ap::ReproducingKernel{:Cubic2D},ξ::NTuple{3,Float64}) = get𝒑₂(ap,ξ[1],ξ[2])
@inline get𝒑₂(::ReproducingKernel{:Cubic2D},ξ::Float64,η::Float64) = (1.,ξ,η)
@inline get∂𝒑₂∂ξ(ap::ReproducingKernel{:Cubic2D},ξ::Any) = (0.,1.,0.)
@inline get∂𝒑₂∂η(ap::ReproducingKernel{:Cubic2D},ξ::Any) = (0.,0.,1.)
@inline get∂²𝒑₂∂ξ²(ap::ReproducingKernel{:Cubic2D},ξ::Any) = (0.,0.,0.)
@inline get∂²𝒑₂∂ξ∂η(ap::ReproducingKernel{:Cubic2D},ξ::Any) = (0.,0.,0.)
@inline get∂²𝒑₂∂η²(ap::ReproducingKernel{:Cubic2D},ξ::Any) = (0.,0.,0.)

## Kernel Function
function get𝜙(ap::ReproducingKernel{𝒑,:□,𝜙},x::Node,Δx::NTuple{3,Float64}) where {𝒑,𝜙}
    rx = abs(Δx[1])/x.s₁
    ry = abs(Δx[2])/x.s₂
    rz = abs(Δx[3])/x.s₃
    wx = get𝜙ᵣ(ap,rx)
    wy = get𝜙ᵣ(ap,ry)
    wz = get𝜙ᵣ(ap,rz)
    return wx*wy*wz
end

function get∇𝜙(ap::ReproducingKernel{𝒑,:□,𝜙},x::Node,Δx::NTuple{3,Float64}) where {𝒑,𝜙}
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

function get∇²𝜙(ap::ReproducingKernel{𝒑,:□,𝜙},x::Node,Δx::NTuple{3,Float64}) where {𝒑,𝜙}
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

function get∇³𝜙(ap::ReproducingKernel{𝒑,:□,𝜙},x::Node,Δx::NTuple{3,Float64}) where {𝒑,𝜙}
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
function get𝜙ᵣ(::ReproducingKernel{𝒑,𝑠,:CubicSpline},r::Float64) where {𝒑,𝑠}
    if r > 1.0
        return 0.0
    elseif r <= 0.5
        return 2/3 - 4*r^2 +  4*r^3
    else
        return 4/3 - 4*r + 4*r^2 - 4*r^3/3
    end
end

function get∂𝜙∂r(::ReproducingKernel{𝒑,𝑠,:CubicSpline},r::Float64) where {𝒑,𝑠}
    if r > 1.0
        return 0.0
    elseif r <= 0.5
        return - 8*r + 12*r^2
    else
        return - 4   + 8*r - 4*r^2
    end
end

function get∂²𝜙∂r²(::ReproducingKernel{𝒑,𝑠,:CubicSpline},r::Float64) where {𝒑,𝑠}
    if r > 1.0
        return 0.0
    elseif r <= 0.5
        return - 8 + 24*r
    else
        return   8 - 8*r
    end
end

function get𝜙ᵣ(::ReproducingKernel{𝒑,𝑠,:QuinticSpline},r::Float64) where {𝒑,𝑠}
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

function get∂𝜙∂r(::ReproducingKernel{𝒑,𝑠,:QuinticSpline},r::Float64) where {𝒑,𝑠}
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

function get∂²𝜙∂r²(::ReproducingKernel{𝒑,𝑠,:QuinticSpline},r::Float64) where {𝒑,𝑠}
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

function get∂³𝜙∂r³(::ReproducingKernel{𝒑,𝑠,:QuinticSpline},r::Float64) where {𝒑,𝑠}
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
function cal𝗠!(ap::ReproducingKernel,x::NTuple{3,Float64})
    𝓒 = ap.𝓒
    𝗠 = ap.𝗠[:∂1]
    n = get𝑛𝒑(ap)
    fill!(𝗠,0.)
    for xᵢ in 𝓒
        Δx = x - xᵢ
        𝒑 = get𝒑(ap,Δx)
        𝜙 = get𝜙(ap,xᵢ,Δx)
        for I in 1:n
            for J in 1:I
                𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
            end
        end
    end
    cholesky!(𝗠)
    inverse!(𝗠)
    UUᵀ!(𝗠)
end

function cal∇₁𝗠!(ap::ReproducingKernel,x::NTuple{3,Float64})
    𝓒 = ap.𝓒
    𝗠 = ap.𝗠[:∂1]
    ∂𝗠∂x = ap.𝗠[:∂x]
    n = get𝑛𝒑(ap)
    fill!(𝗠,0.)
    fill!(∂𝗠∂x,0.)
    for xᵢ in 𝓒
        Δx = x - xᵢ
        𝒑, ∂𝒑∂x = get∇₁𝒑(ap,Δx)
        𝜙, ∂𝜙∂x = get∇₁𝜙(ap,xᵢ,Δx)
        for I in 1:n
            for J in 1:I
                𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
                ∂𝗠∂x[I,J] += ∂𝜙∂x*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂x[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂x[J]
            end
        end
    end
    cholesky!(𝗠)
    U = inverse!(𝗠)
    ∂𝗠⁻¹∂x = - UUᵀAUUᵀ!(∂𝗠∂x,U)
    𝗠⁻¹ = UUᵀ!(U)
    return 𝗠⁻¹, ∂𝗠⁻¹∂x
end

function cal∇₂𝗠!(ap::ReproducingKernel,x::NTuple{3,Float64})
    𝓒 = ap.𝓒
    𝗠 = ap.𝗠[:∂1]
    ∂𝗠∂x = ap.𝗠[:∂x]
    ∂𝗠∂y = ap.𝗠[:∂y]
    n = get𝑛𝒑(ap)
    fill!(𝗠,0.)
    fill!(∂𝗠∂x,0.)
    fill!(∂𝗠∂y,0.)
    for xᵢ in 𝓒
        Δx = x - xᵢ
        𝒑, ∂𝒑∂x, ∂𝒑∂y = get∇₂𝒑(ap,Δx)
        𝜙, ∂𝜙∂x, ∂𝜙∂y = get∇₂𝜙(ap,xᵢ,Δx)
        for I in 1:n
            for J in 1:I
                𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
                ∂𝗠∂x[I,J] += ∂𝜙∂x*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂x[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂x[J]
                ∂𝗠∂y[I,J] += ∂𝜙∂y*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂y[J]
            end
        end
    end
    cholesky!(𝗠)
    U = inverse!(𝗠)
    ∂𝗠⁻¹∂x = - UUᵀAUUᵀ!(∂𝗠∂x,U)
    ∂𝗠⁻¹∂y = - UUᵀAUUᵀ!(∂𝗠∂y,U)
    𝗠⁻¹ = UUᵀ!(U)
    return 𝗠⁻¹, ∂𝗠⁻¹∂x, ∂𝗠⁻¹∂y
end

function cal∇₃𝗠!(ap::ReproducingKernel,x::NTuple{3,Float64})
    𝓒 = ap.𝓒
    𝗠 = ap.𝗠[:∂1]
    ∂𝗠∂x = ap.𝗠[:∂x]
    ∂𝗠∂y = ap.𝗠[:∂y]
    ∂𝗠∂z = ap.𝗠[:∂z]
    n = get𝑛𝒑(ap)
    fill!(𝗠,0.)
    fill!(∂𝗠∂x,0.)
    fill!(∂𝗠∂y,0.)
    fill!(∂𝗠∂z,0.)
    for xᵢ in 𝓒
        Δx = x - xᵢ
        𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂𝒑∂z = get∇₃𝒑(ap,Δx)
        𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂𝜙∂z = get∇₃𝜙(ap,xᵢ,Δx)
        for I in 1:n
            for J in 1:I
                𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
                ∂𝗠∂x[I,J] += ∂𝜙∂x*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂x[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂x[J]
                ∂𝗠∂y[I,J] += ∂𝜙∂y*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂y[J]
                ∂𝗠∂z[I,J] += ∂𝜙∂z*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂z[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂z[J]
            end
        end
    end
    cholesky!(𝗠)
    U = inverse!(𝗠)
    ∂𝗠⁻¹∂x = - UUᵀAUUᵀ!(∂𝗠∂x,U)
    ∂𝗠⁻¹∂y = - UUᵀAUUᵀ!(∂𝗠∂y,U)
    ∂𝗠⁻¹∂z = - UUᵀAUUᵀ!(∂𝗠∂z,U)
    𝗠⁻¹ = UUᵀ!(U)
    return 𝗠⁻¹, ∂𝗠⁻¹∂x, ∂𝗠⁻¹∂y, ∂𝗠⁻¹∂z
end

function cal∇²₁𝗠!(ap::ReproducingKernel,x::NTuple{3,Float64})
    𝓒 = ap.𝓒
    𝗠 = ap.𝗠[:∂1]
    ∂𝗠∂x = ap.𝗠[:∂x]
    ∂²𝗠∂x² = ap.𝗠[:∂x²]
    n = get𝑛𝒑(ap)
    fill!(𝗠,0.)
    fill!(∂𝗠∂x,0.)
    fill!(∂²𝗠∂x²,0.)
    for xᵢ in 𝓒
        Δx = x - xᵢ
        𝒑, ∂𝒑∂x, ∂²𝒑∂x² = get∇²₁𝒑(ap,Δx)
        𝜙, ∂𝜙∂x, ∂²𝜙∂x² = get∇²₁𝜙(ap,xᵢ,Δx)
        for I in 1:n
            for J in I:n
                𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
                ∂𝗠∂x[I,J] += ∂𝜙∂x*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂x[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂x[J]
                ∂²𝗠∂x²[I,J] += ∂²𝜙∂x²*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂x²[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂x²[J] + 2.0*∂𝜙∂x*∂𝒑∂x[I]*𝒑[J] + 2.0*∂𝜙∂x*𝒑[I]*∂𝒑∂x[J] + 2.0*𝜙*∂𝒑∂x[I]*∂𝒑∂x[J]
            end
        end
    end
    cholesky!(𝗠)
    U = inverse!(𝗠)
    Uᵀ∂𝗠∂xU = UᵀAU!(∂𝗠∂x,U)
    Uᵀ∂²𝗠∂x²U = UᵀAU!(∂²𝗠∂x²,U)
    for i in 1:n
        for j in 1:i
            for k in 1:n
                Uᵀ∂²𝗠∂x²U[i,j] -= 2*Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂𝗠∂xU[k,j]
            end
        end
    end

    ∂²𝗠⁻¹∂x² = - UAUᵀ!(Uᵀ∂²𝗠∂x²U,U)
    ∂𝗠⁻¹∂x = - UAUᵀ!(Uᵀ∂𝗠∂xU,U)
    𝗠⁻¹ = UUᵀ!(U)
    return 𝗠⁻¹, ∂𝗠⁻¹∂x, ∂²𝗠⁻¹∂x²
end

function cal∇²₂𝗠!(ap::ReproducingKernel,x::NTuple{3,Float64})
    𝓒 = ap.𝓒
    𝗠 = ap.𝗠[:∂1]
    ∂𝗠∂x = ap.𝗠[:∂x]
    ∂𝗠∂y = ap.𝗠[:∂y]
    ∂²𝗠∂x² = ap.𝗠[:∂x²]
    ∂²𝗠∂y² = ap.𝗠[:∂y²]
    ∂²𝗠∂x∂y = ap.𝗠[:∂x∂y]
    n = get𝑛𝒑(ap)
    fill!(𝗠,0.)
    fill!(∂𝗠∂x,0.)
    fill!(∂𝗠∂y,0.)
    fill!(∂²𝗠∂x²,0.)
    fill!(∂²𝗠∂y²,0.)
    fill!(∂²𝗠∂x∂y,0.)
    for xᵢ in 𝓒
        Δx = x - xᵢ
        𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂²𝒑∂x², ∂²𝒑∂x∂y, ∂²𝒑∂y² = get∇²₂𝒑(ap,Δx)
        𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂²𝜙∂x², ∂²𝜙∂x∂y, ∂²𝜙∂y² = get∇²₂𝜙(ap,xᵢ,Δx)
        for I in 1:n
            for J in 1:I
                𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
                ∂𝗠∂x[I,J] += ∂𝜙∂x*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂x[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂x[J]
                ∂𝗠∂y[I,J] += ∂𝜙∂y*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂y[J]
                ∂²𝗠∂x²[I,J] += ∂²𝜙∂x²*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂x²[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂x²[J] + 2.0*∂𝜙∂x*∂𝒑∂x[I]*𝒑[J] + 2.0*∂𝜙∂x*𝒑[I]*∂𝒑∂x[J] + 2.0*𝜙*∂𝒑∂x[I]*∂𝒑∂x[J]
                ∂²𝗠∂y²[I,J] += ∂²𝜙∂y²*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂y²[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂y²[J] + 2.0*∂𝜙∂y*∂𝒑∂y[I]*𝒑[J] + 2.0*∂𝜙∂y*𝒑[I]*∂𝒑∂y[J] + 2.0*𝜙*∂𝒑∂y[I]*∂𝒑∂y[J]
                ∂²𝗠∂x∂y[I,J] += ∂²𝜙∂x∂y*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂x∂y[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂x∂y[J] + ∂𝜙∂x*∂𝒑∂y[I]*𝒑[J] + ∂𝜙∂y*∂𝒑∂x[I]*𝒑[J] + ∂𝜙∂x*𝒑[I]*∂𝒑∂y[J] + ∂𝜙∂y*𝒑[I]*∂𝒑∂x[J] + 𝜙*∂𝒑∂x[I]*∂𝒑∂y[J] + 𝜙*∂𝒑∂y[I]*∂𝒑∂x[J]
            end
        end
    end
    cholesky!(𝗠)
    U = inverse!(𝗠)
    Uᵀ∂𝗠∂xU = UᵀAU!(∂𝗠∂x,U)
    Uᵀ∂𝗠∂yU = UᵀAU!(∂𝗠∂y,U)
    Uᵀ∂²𝗠∂x²U = UᵀAU!(∂²𝗠∂x²,U)
    Uᵀ∂²𝗠∂y²U = UᵀAU!(∂²𝗠∂y²,U)
    Uᵀ∂²𝗠∂x∂yU = UᵀAU!(∂²𝗠∂x∂y,U)
    for i in 1:n
        for j in 1:i
            for k in 1:n
                Uᵀ∂²𝗠∂x²U[i,j] -= 2*Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂𝗠∂xU[k,j]
                Uᵀ∂²𝗠∂y²U[i,j] -= 2*Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂𝗠∂yU[k,j]
                Uᵀ∂²𝗠∂x∂yU[i,j] -= Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂𝗠∂yU[k,j] + Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂𝗠∂xU[k,j]
            end
        end
    end

    ∂²𝗠⁻¹∂x² = - UAUᵀ!(Uᵀ∂²𝗠∂x²U,U)
    ∂²𝗠⁻¹∂y² = - UAUᵀ!(Uᵀ∂²𝗠∂y²U,U)
    ∂²𝗠⁻¹∂x∂y = - UAUᵀ!(Uᵀ∂²𝗠∂x∂yU,U)
    ∂𝗠⁻¹∂x = - UAUᵀ!(Uᵀ∂𝗠∂xU,U)
    ∂𝗠⁻¹∂y = - UAUᵀ!(Uᵀ∂𝗠∂yU,U)
    𝗠⁻¹ = UUᵀ!(U)
    return 𝗠⁻¹, ∂𝗠⁻¹∂x, ∂𝗠⁻¹∂y, ∂²𝗠⁻¹∂x², ∂²𝗠⁻¹∂x∂y, ∂²𝗠⁻¹∂y²
end

function cal∇²₃𝗠!(ap::ReproducingKernel,x::NTuple{3,Float64})
    𝓒 = ap.𝓒
    𝗠 = ap.𝗠[:∂1]
    ∂𝗠∂x = ap.𝗠[:∂x]
    ∂𝗠∂y = ap.𝗠[:∂y]
    ∂𝗠∂z = ap.𝗠[:∂z]
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
        𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂²𝒑∂x², ∂²𝒑∂x∂y, ∂²𝒑∂y², ∂𝒑∂z, ∂²𝒑∂x∂z, ∂²𝒑∂y∂z, ∂²𝒑∂z² = get∇²₃𝒑(ap,Δx)
        𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂²𝜙∂x², ∂²𝜙∂x∂y, ∂²𝜙∂y², ∂𝜙∂z, ∂²𝜙∂x∂z, ∂²𝜙∂y∂z, ∂²𝜙∂z² = get∇²₃𝜙(ap,xᵢ,Δx)
        for I in 1:n
            for J in 1:I
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
    U = inverse!(𝗠)
    Uᵀ∂𝗠∂xU = UᵀAU!(∂𝗠∂x,U)
    Uᵀ∂𝗠∂yU = UᵀAU!(∂𝗠∂y,U)
    Uᵀ∂𝗠∂zU = UᵀAU!(∂𝗠∂z,U)
    Uᵀ∂²𝗠∂x²U = UᵀAU!(∂²𝗠∂x²,U)
    Uᵀ∂²𝗠∂y²U = UᵀAU!(∂²𝗠∂y²,U)
    Uᵀ∂²𝗠∂z²U = UᵀAU!(∂²𝗠∂z²,U)
    Uᵀ∂²𝗠∂x∂yU = UᵀAU!(∂²𝗠∂x∂y,U)
    Uᵀ∂²𝗠∂x∂zU = UᵀAU!(∂²𝗠∂x∂z,U)
    Uᵀ∂²𝗠∂y∂zU = UᵀAU!(∂²𝗠∂y∂z,U)
    for i in 1:n
        for j in 1:i
            for k in 1:n
                Uᵀ∂²𝗠∂x²U[i,j] -= 2*Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂𝗠∂xU[k,j]
                Uᵀ∂²𝗠∂y²U[i,j] -= 2*Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂𝗠∂yU[k,j]
                Uᵀ∂²𝗠∂z²U[i,j] -= 2*Uᵀ∂𝗠∂zU[i,k]*Uᵀ∂𝗠∂zU[k,j]
                Uᵀ∂²𝗠∂x∂yU[i,j] -= Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂𝗠∂yU[k,j] + Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂𝗠∂xU[k,j]
                Uᵀ∂²𝗠∂x∂zU[i,j] -= Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂𝗠∂zU[k,j] + Uᵀ∂𝗠∂zU[i,k]*Uᵀ∂𝗠∂xU[k,j]
                Uᵀ∂²𝗠∂y∂zU[i,j] -= Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂𝗠∂zU[k,j] + Uᵀ∂𝗠∂zU[i,k]*Uᵀ∂𝗠∂yU[k,j]
            end
        end
    end

    ∂²𝗠⁻¹∂x² = - UAUᵀ!(Uᵀ∂²𝗠∂x²U,U)
    ∂²𝗠⁻¹∂y² = - UAUᵀ!(Uᵀ∂²𝗠∂y²U,U)
    ∂²𝗠⁻¹∂z² = - UAUᵀ!(Uᵀ∂²𝗠∂z²U,U)
    ∂²𝗠⁻¹∂x∂y = - UAUᵀ!(Uᵀ∂²𝗠∂x∂yU,U)
    ∂²𝗠⁻¹∂x∂z = - UAUᵀ!(Uᵀ∂²𝗠∂x∂zU,U)
    ∂²𝗠⁻¹∂y∂z = - UAUᵀ!(Uᵀ∂²𝗠∂y∂zU,U)
    ∂𝗠⁻¹∂x = - UAUᵀ!(Uᵀ∂𝗠∂xU,U)
    ∂𝗠⁻¹∂y = - UAUᵀ!(Uᵀ∂𝗠∂yU,U)
    ∂𝗠⁻¹∂z = - UAUᵀ!(Uᵀ∂𝗠∂zU,U)
    𝗠⁻¹ = UUᵀ!(U)
    return 𝗠⁻¹, ∂𝗠⁻¹∂x, ∂𝗠⁻¹∂y, ∂²𝗠⁻¹∂x², ∂²𝗠⁻¹∂x∂y, ∂²𝗠⁻¹∂y², ∂𝗠⁻¹∂z, ∂²𝗠⁻¹∂x∂z, ∂²𝗠⁻¹∂y∂z, ∂²𝗠⁻¹∂z²
end

function cal∇³₁𝗠!(ap::ReproducingKernel,x::NTuple{3,Float64})
    𝓒 = ap.𝓒
    𝗠 = ap.𝗠[:∂1]
    ∂𝗠∂x = ap.𝗠[:∂x]
    ∂𝗠∂y = ap.𝗠[:∂y]
    ∂²𝗠∂x² = ap.𝗠[:∂x²]
    ∂²𝗠∂x∂y = ap.𝗠[:∂x∂y]
    ∂²𝗠∂y² = ap.𝗠[:∂y²]
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
        𝒑, ∂𝒑∂x, ∂²𝒑∂x², ∂³𝒑∂x³ = get∇³₁𝒑(ap,Δx)
        𝜙, ∂𝜙∂x, ∂²𝜙∂x², ∂³𝜙∂x³ = get∇³₁𝜙(ap,xᵢ,Δx)
        for I in 1:n
            for J in 1:I
                𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
                ∂𝗠∂x[I,J] += ∂𝜙∂x*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂x[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂x[J]
                ∂²𝗠∂x²[I,J] += ∂²𝜙∂x²*𝒑[I]*𝒑[J] + 𝜙*∂²𝒑∂x²[I]*𝒑[J] + 𝜙*𝒑[I]*∂²𝒑∂x²[J] + 2.0*∂𝜙∂x*∂𝒑∂x[I]*𝒑[J] + 2.0*∂𝜙∂x*𝒑[I]*∂𝒑∂x[J] + 2.0*𝜙*∂𝒑∂x[I]*∂𝒑∂x[J]
                ∂³𝗠∂x³[I,J] += ∂³𝜙∂x³*𝒑[I]*𝒑[J] + 𝜙*∂³𝒑∂x³[I]*𝒑[J] + 𝜙*𝒑[I]*∂³𝒑∂x³[J] + 3.0*∂²𝜙∂x²*∂𝒑∂x[I]*𝒑[J] + 3.0*∂𝜙∂x*∂²𝒑∂x²[I]*𝒑[J] + 3.0*∂²𝜙∂x²*𝒑[I]*∂𝒑∂x[J] + 3.0*∂𝜙∂x*𝒑[I]*∂²𝒑∂x²[J] + 3.0*𝜙*∂²𝒑∂x²[I]*∂𝒑∂x[J] + 3.0*𝜙*∂𝒑∂x[I]*∂²𝒑∂x²[J] + 6.0*∂𝜙∂x*∂𝒑∂x[I]*∂𝒑∂x[J]
            end
        end
    end
    cholesky!(𝗠)
    U = inverse!(𝗠)
    Uᵀ∂𝗠∂xU = UᵀAU!(∂𝗠∂x,U)
    Uᵀ∂²𝗠∂x²U = UᵀAU!(∂²𝗠∂x²,U)
    Uᵀ∂³𝗠∂x³U = UᵀAU!(∂³𝗠∂x³,U)

    for i in 1:n
        for j in 1:i
            for k in 1:n
                Uᵀ∂³𝗠∂x³U[i,j] -= 3*Uᵀ∂²𝗠∂x²U[i,k]*Uᵀ∂𝗠∂xU[k,j]
            end
        end
        for j in 1:i
            for k in 1:n
                Uᵀ∂²𝗠∂x²U[i,j] -= 2*Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂𝗠∂xU[k,j]
            end
        end
    end
    for i in 1:n
        for j in 1:i
            for k in 1:n
                Uᵀ∂³𝗠∂x³U[i,j] -= 3*Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂²𝗠∂x²U[k,j]
            end
        end
    end

    ∂³𝗠⁻¹∂x³ = - UAUᵀ!(Uᵀ∂³𝗠∂x³U,U)
    ∂²𝗠⁻¹∂x² = - UAUᵀ!(Uᵀ∂²𝗠∂x²U,U)
    ∂𝗠⁻¹∂x = - UAUᵀ!(Uᵀ∂𝗠∂xU,U)
    𝗠⁻¹ = UUᵀ!(U)

    return 𝗠⁻¹, ∂𝗠⁻¹∂x, ∂²𝗠⁻¹∂x², ∂³𝗠⁻¹∂x³
end

function cal∇³₂𝗠!(ap::ReproducingKernel,x::NTuple{3,Float64})
    𝓒 = ap.𝓒
    𝗠 = ap.𝗠[:∂1]
    ∂𝗠∂x = ap.𝗠[:∂x]
    ∂𝗠∂y = ap.𝗠[:∂y]
    ∂²𝗠∂x² = ap.𝗠[:∂x²]
    ∂²𝗠∂x∂y = ap.𝗠[:∂x∂y]
    ∂²𝗠∂y² = ap.𝗠[:∂y²]
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
        𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂²𝒑∂x², ∂²𝒑∂x∂y, ∂²𝒑∂y², ∂³𝒑∂x³, ∂³𝒑∂x²∂y, ∂³𝒑∂x∂y², ∂³𝒑∂y³ = get∇³₂𝒑(ap,Δx)
        𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂²𝜙∂x², ∂²𝜙∂x∂y, ∂²𝜙∂y², ∂³𝜙∂x³, ∂³𝜙∂x²∂y, ∂³𝜙∂x∂y², ∂³𝜙∂y³ = get∇³₂𝜙(ap,xᵢ,Δx)
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
    U = inverse!(𝗠)
    Uᵀ∂𝗠∂xU = UᵀAU!(∂𝗠∂x,U)
    Uᵀ∂𝗠∂yU = UᵀAU!(∂𝗠∂y,U)
    Uᵀ∂²𝗠∂x²U = UᵀAU!(∂²𝗠∂x²,U)
    Uᵀ∂²𝗠∂y²U = UᵀAU!(∂²𝗠∂y²,U)
    Uᵀ∂²𝗠∂x∂yU = UᵀAU!(∂²𝗠∂x∂y,U)
    Uᵀ∂³𝗠∂x³U = UᵀAU!(∂³𝗠∂x³,U)
    Uᵀ∂³𝗠∂x²∂yU = UᵀAU!(∂³𝗠∂x²∂y,U)
    Uᵀ∂³𝗠∂x∂y²U = UᵀAU!(∂³𝗠∂x∂y²,U)
    Uᵀ∂³𝗠∂y³U = UᵀAU!(∂³𝗠∂y³,U)

    for i in 1:n
        for j in 1:i
            for k in 1:n
                Uᵀ∂³𝗠∂x³U[i,j] -= 3*Uᵀ∂²𝗠∂x²U[i,k]*Uᵀ∂𝗠∂xU[k,j]
                Uᵀ∂³𝗠∂x²∂yU[i,j] -= 2*Uᵀ∂²𝗠∂x∂yU[i,k]*Uᵀ∂𝗠∂xU[k,j]+Uᵀ∂²𝗠∂x²U[i,k]*Uᵀ∂𝗠∂yU[k,j]
                Uᵀ∂³𝗠∂x∂y²U[i,j] -= 2*Uᵀ∂²𝗠∂x∂yU[i,k]*Uᵀ∂𝗠∂yU[k,j]+Uᵀ∂²𝗠∂y²U[i,k]*Uᵀ∂𝗠∂xU[k,j]
                Uᵀ∂³𝗠∂y³U[i,j] -= 3*Uᵀ∂²𝗠∂y²U[i,k]*Uᵀ∂𝗠∂yU[k,j]
            end
        end
        for j in 1:i
            for k in 1:n
                Uᵀ∂²𝗠∂x²U[i,j] -= 2*Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂𝗠∂xU[k,j]
                Uᵀ∂²𝗠∂y²U[i,j] -= 2*Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂𝗠∂yU[k,j]
                Uᵀ∂²𝗠∂x∂yU[i,j] -= Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂𝗠∂yU[k,j] + Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂𝗠∂xU[k,j]
            end
        end
    end
    for i in 1:n
        for j in 1:i
            for k in 1:n
                Uᵀ∂³𝗠∂x³U[i,j] -= 3*Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂²𝗠∂x²U[k,j]
                Uᵀ∂³𝗠∂x²∂yU[i,j] -= 2*Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂²𝗠∂x∂yU[k,j]+Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂²𝗠∂x²U[k,j]
                Uᵀ∂³𝗠∂x∂y²U[i,j] -= 2*Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂²𝗠∂x∂yU[k,j]+Uᵀ∂𝗠∂xU[i,k]*Uᵀ∂²𝗠∂y²U[k,j]
                Uᵀ∂³𝗠∂y³U[i,j] -= 3*Uᵀ∂𝗠∂yU[i,k]*Uᵀ∂²𝗠∂y²U[k,j]
            end
        end
    end

    ∂³𝗠⁻¹∂x³ = - UAUᵀ!(Uᵀ∂³𝗠∂x³U,U)
    ∂³𝗠⁻¹∂x²∂y = - UAUᵀ!(Uᵀ∂³𝗠∂x²∂yU,U)
    ∂³𝗠⁻¹∂x∂y² = - UAUᵀ!(Uᵀ∂³𝗠∂x∂y²U,U)
    ∂³𝗠⁻¹∂y³ = - UAUᵀ!(Uᵀ∂³𝗠∂y³U,U)
    ∂²𝗠⁻¹∂x² = - UAUᵀ!(Uᵀ∂²𝗠∂x²U,U)
    ∂²𝗠⁻¹∂y² = - UAUᵀ!(Uᵀ∂²𝗠∂y²U,U)
    ∂²𝗠⁻¹∂x∂y = - UAUᵀ!(Uᵀ∂²𝗠∂x∂yU,U)
    ∂𝗠⁻¹∂x = - UAUᵀ!(Uᵀ∂𝗠∂xU,U)
    ∂𝗠⁻¹∂y = - UAUᵀ!(Uᵀ∂𝗠∂yU,U)
    𝗠⁻¹ = UUᵀ!(U)

    return 𝗠⁻¹, ∂𝗠⁻¹∂x, ∂𝗠⁻¹∂y, ∂²𝗠⁻¹∂x², ∂²𝗠⁻¹∂x∂y, ∂²𝗠⁻¹∂y², ∂³𝗠⁻¹∂x³, ∂³𝗠⁻¹∂x²∂y, ∂³𝗠⁻¹∂x∂y², ∂³𝗠⁻¹∂y³
end
