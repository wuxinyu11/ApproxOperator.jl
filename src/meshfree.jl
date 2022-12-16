"""
SymMat
"""
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

"""
Spatial Partition
RegularGrid 
"""
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

function RegularGrid(nodes::Vector{T};n::Int=1,γ::Int=1) where T<:AbstractNode
    node = nodes[1]
    x = getfield(node,:data)[:x][2]
    y = getfield(node,:data)[:y][2]
    z = getfield(node,:data)[:z][2]
    return RegularGrid(x,y,z,n=n,γ=γ)
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
        𝓒 = ap.𝓒; 𝓖 = ap.𝓖
        indices = Set{Int}()
        for 𝒙 in 𝓒
            union!(indices,sp(𝒙.x,𝒙.y,𝒙.z))
        end
        union!(𝓒,(Node(i,getfield(𝓒[1],:data)) for i in indices))
    end
    function (sp::t)(aps::Vector{T}) where T<:AbstractElement
        for ap in aps
            sp(ap)
        end
    end
    function (sp::t)(ap::T,nodes::Vector{Node}) where T<:AbstractElement
        𝓒 = ap.𝓒; 𝓖 = ap.𝓖
        indices = Set{Int}()
        for 𝒙 in 𝓒
            union!(indices,sp(𝒙.x,𝒙.y,𝒙.z))
        end
        return [nodes[i] for i in indices]
    end
    function (sp::t)(apss::Any...)
        for aps in apss
            sp(aps)
        end
    end
end

"""
ReproducingKernel
"""
struct ReproducingKernel{𝑝,𝑠,𝜙,T}<:AbstractElement{T}
    𝓒::Vector{Node}
    𝓖::Vector{SNode}
end

function get𝗠(ap::ReproducingKernel,s::Symbol)
    n = get𝑛𝒑(ap)
    data = getfield(ap.𝓖[1],:data)
    fill!(data[s][2],0.)
    return SymMat(n,data[s][2])
end
function get𝗚(ap::ReproducingKernel,s::Symbol)
    n = get𝑛𝒑₁(ap)
    data = getfield(ap.𝓖[1],:data)
    fill!(data[s][2],0.)
    return SymMat(n,data[s][2])
end

"""
Basis function
"""
## Basis Function
@inline get∇₁𝒑(ap::ReproducingKernel,x::Any) = get𝒑(ap,x), get∂𝒑∂x(ap,x)
@inline get∇₂𝒑(ap::ReproducingKernel,x::Any) = get𝒑(ap,x), get∂𝒑∂x(ap,x), get∂𝒑∂y(ap,x)
@inline get∇𝒑(ap::ReproducingKernel,x::Any) = get𝒑(ap,x), get∂𝒑∂x(ap,x), get∂𝒑∂y(ap,x), get∂𝒑∂z(ap,x)
@inline get∇²₁𝒑(ap::ReproducingKernel,x::Any) = get𝒑(ap,x), get∂𝒑∂x(ap,x), get∂²𝒑∂x²(ap,x)
@inline get∇²₂𝒑(ap::ReproducingKernel,x::Any) = get𝒑(ap,x), get∂𝒑∂x(ap,x), get∂𝒑∂y(ap,x), get∂²𝒑∂x²(ap,x), get∂²𝒑∂x∂y(ap,x), get∂²𝒑∂y²(ap,x)
@inline get∇̃²₂𝒑(ap::ReproducingKernel,x::Any) = get∂²𝒑∂x²(ap,x), get∂²𝒑∂x∂y(ap,x), get∂²𝒑∂x∂y(ap,x), get∂²𝒑∂y²(ap,x)
@inline get∇²𝒑(ap::ReproducingKernel,x::Any) = get𝒑(ap,x), get∂𝒑∂x(ap,x), get∂𝒑∂y(ap,x), get∂²𝒑∂x²(ap,x), get∂²𝒑∂x∂y(ap,x), get∂²𝒑∂y²(ap,x), get∂𝒑∂z(ap,x), get∂²𝒑∂x∂z(ap,x), get∂²𝒑∂y∂z(ap,x), get∂²𝒑∂z²(ap,x)
@inline get∇³₁𝒑(ap::ReproducingKernel,x::Any) = get𝒑(ap,x), get∂𝒑∂x(ap,x), get∂²𝒑∂x²(ap,x), get∂³𝒑∂x³(ap,x)
@inline get∇³𝒑(ap::ReproducingKernel,x::Any) = get𝒑(ap,x), get∂𝒑∂x(ap,x), get∂𝒑∂y(ap,x), get∂²𝒑∂x²(ap,x), get∂²𝒑∂x∂y(ap,x), get∂²𝒑∂y²(ap,x), get∂³𝒑∂x³(ap,x), get∂³𝒑∂x²∂y(ap,x), get∂³𝒑∂x∂y²(ap,x), get∂³𝒑∂y³(ap,x)
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
@inline get𝒑₁(ap::ReproducingKernel{:Quadratic1D},ξ::SNode) = get𝒑₁(ap,ξ.ξ)
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
@inline get𝒑₁(ap::ReproducingKernel{:Cubic1D},ξ::SNode) = get𝒑₁(ap,ξ.ξ)
@inline get𝒑₁(::ReproducingKernel{:Cubic1D},ξ::Float64) = (1.0,0.5*(1.0-ξ),0.25*(1.0-ξ)^2)
@inline get∂𝒑₁∂ξ(ap::ReproducingKernel{:Cubic1D},ξ::SNode) = get∂𝒑₁∂ξ(ap,ξ.ξ)
@inline get∂𝒑₁∂ξ(::ReproducingKernel{:Cubic1D},ξ::Float64) = (0.,1.0,(1.0-ξ))

# ------------ Linear2D ---------------
@inline get𝑛𝒑(::ReproducingKernel{:Linear2D}) = 3
@inline get𝒑(::ReproducingKernel{:Linear2D},x::NTuple{3,Float64}) = (1.,x[1],x[2])
@inline get∂𝒑∂x(::ReproducingKernel{:Linear2D}, ::Any) = (0.,1.,0.)
@inline get∂𝒑∂y(::ReproducingKernel{:Linear2D}, ::Any) = (0.,0.,1.)
@inline get∂𝒑∂z(::ReproducingKernel{:Linear2D}, ::Any) = (0.,0.,0.)

@inline get𝑛𝒑₁(::ReproducingKernel{:Linear2D}) = 1
@inline get𝒑₁(ap::ReproducingKernel{:Linear2D},ξ::SNode) = get𝒑₁(ap,ξ.ξ,ξ.η)
@inline get𝒑₁(::ReproducingKernel{:Linear2D},::Any,::Any) = (1.,)
@inline get∂𝒑₁∂ξ(ap::ReproducingKernel{:Linear2D},ξ::SNode) = get∂𝒑₁∂ξ(ap,ξ.ξ,ξ.η)
@inline get∂𝒑₁∂ξ(::ReproducingKernel{:Linear2D},::Any,::Any) = (0.,)
@inline get∂𝒑₁∂η(ap::ReproducingKernel{:Linear2D},ξ::SNode) = get∂𝒑₁∂η(ap,ξ.ξ,ξ.η)
@inline get∂𝒑₁∂η(::ReproducingKernel{:Linear2D},::Any,::Any) = (0.,)

@inline get𝑛𝒑₂(::ReproducingKernel{:Linear2D}) = 0
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
@inline get𝒑₁(ap::ReproducingKernel{:Quadratic2D},ξ::SNode) = get𝒑₁(ap,ξ.ξ,ξ.η)
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
@inline get𝒑₁(ap::ReproducingKernel{:Cubic2D},ξ::SNode) = get𝒑₁(ap,ξ.ξ,ξ.η)
@inline get𝒑₁(::ReproducingKernel{:Cubic2D},ξ::Float64,η::Float64) = (1.,ξ,η,ξ^2,ξ*η,η^2)
@inline get∂𝒑₁∂ξ(ap::ReproducingKernel{:Cubic2D},ξ::SNode) = get∂𝒑₁∂ξ(ap,ξ.ξ,ξ.η)
@inline get∂𝒑₁∂ξ(::ReproducingKernel{:Cubic2D},ξ::Float64,η::Float64) = (0.,1.,0.,2.0*ξ,η,0.)
@inline get∂𝒑₁∂η(ap::ReproducingKernel{:Cubic2D},ξ::SNode) = get∂𝒑₁∂η(ap,ξ.ξ,ξ.η)
@inline get∂𝒑₁∂η(::ReproducingKernel{:Cubic2D},ξ::Float64,η::Float64) = (0.,0.,1.,0.,ξ,2.0*η)

@inline get𝑛𝒑₂(::ReproducingKernel{:Cubic2D}) = 3
@inline get𝒑₂(ap::ReproducingKernel{:Cubic2D},ξ::SNode) = get𝒑₂(ap,ξ.ξ,ξ.η)
@inline get𝒑₂(ap::ReproducingKernel{:Cubic2D},ξ::NTuple{3,Float64}) = get𝒑₂(ap,ξ[1],ξ[2])
@inline get𝒑₂(::ReproducingKernel{:Cubic2D},ξ::Float64,η::Float64) = (1.,ξ,η)
@inline get∂𝒑₂∂ξ(ap::ReproducingKernel{:Cubic2D},ξ::Any) = (0.,1.,0.)
@inline get∂𝒑₂∂η(ap::ReproducingKernel{:Cubic2D},ξ::Any) = (0.,0.,1.)
@inline get∂²𝒑₂∂ξ²(ap::ReproducingKernel{:Cubic2D},ξ::Any) = (0.,0.,0.)
@inline get∂²𝒑₂∂ξ∂η(ap::ReproducingKernel{:Cubic2D},ξ::Any) = (0.,0.,0.)
@inline get∂²𝒑₂∂η²(ap::ReproducingKernel{:Cubic2D},ξ::Any) = (0.,0.,0.)

# ------------ Quartic2D ---------------
@inline get𝑛𝒑(::ReproducingKernel{:Quartic2D}) = 15
@inline get𝒑(::ReproducingKernel{:Quartic2D},x::NTuple{3,Float64}) =
(
    1., x[1], x[2], x[1]^2, x[1]*x[2], x[2]^2, x[1]^3, x[1]^2*x[2], x[1]*x[2]^2, x[2]^3, x[1]^4, x[1]^3*x[2], x[1]^2*x[2]^2, x[1]*x[2]^3, x[2]^4
)
@inline get∂𝒑∂x(::ReproducingKernel{:Quartic2D},x::NTuple{3,Float64}) =
(
    0., 1., 0., 2*x[1], x[2], 0., 3*x[1]^2, 2*x[1]*x[2], x[2]^2, 0., 4.0*x[1]^3, 3.0*x[1]^2*x[2], 2.0*x[1]*x[2]^2, x[2]^3, 0.
)
@inline get∂𝒑∂y(::ReproducingKernel{:Quartic2D},x::NTuple{3,Float64}) =
(
    0., 0., 1., 0., x[1], 2*x[2], 0., x[1]^2, 2*x[1]*x[2], 3*x[2]^2, 0.0, x[1]^3, 2.0*x[1]^2*x[2], 3.0*x[1]*x[2]^2, 4.0*x[2]^3
)
@inline get∂²𝒑∂x²(::ReproducingKernel{:Quartic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 2., 0., 0., 6*x[1], 2*x[2], 0., 0., 12.0*x[1]^2, 6.0*x[1]*x[2], 2.0*x[2]^2, 0.0, 0.0
)
@inline get∂²𝒑∂x∂y(::ReproducingKernel{:Quartic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 1., 0., 0., 2*x[1], 2*x[2], 0., 0.0, 3.0*x[1]^2, 4.0*x[1]*x[2], 3.0*x[2]^2, 0.0
)
@inline get∂²𝒑∂y²(::ReproducingKernel{:Quartic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 2., 0., 0., 2*x[1], 6*x[2], 0.0, 0.0, 2.0*x[1]^2, 6.0*x[1]*x[2], 12.0*x[2]^2
)
@inline get∂𝒑∂z(::ReproducingKernel{:Quartic2D},::Any) =
(
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
)
@inline get∂²𝒑∂x∂z(::ReproducingKernel{:Quartic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
)
@inline get∂²𝒑∂y∂z(::ReproducingKernel{:Quartic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
)
@inline get∂²𝒑∂z²(::ReproducingKernel{:Quartic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
)
@inline get∂³𝒑∂x³(::ReproducingKernel{:Quartic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 0., 6., 0., 0., 0., 24.0*x[1], 6.0*x[2], 0.0, 0.0, 0.0
)
@inline get∂³𝒑∂x²∂y(::ReproducingKernel{:Quartic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0., 6.0*x[1], 4.0*x[2], 0., 0.
)
@inline get∂³𝒑∂x∂y²(::ReproducingKernel{:Quartic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 0., 0., 0., 2., 0., 0., 0., 4.0*x[1], 6.0*x[2],0.
)
@inline get∂³𝒑∂y³(::ReproducingKernel{:Quartic2D},x::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 0., 0., 0., 0., 6., 0., 0., 0., 6.0*x[1], 24.0*x[2]
)

@inline get𝑛𝒑₁(::ReproducingKernel{:Quartic2D}) = 10
@inline get𝒑₁(ap::ReproducingKernel{:Quartic2D},ξ::SNode) = get𝒑₁(ap,ξ.ξ,ξ.η)
@inline get𝒑₁(::ReproducingKernel{:Quartic2D},ξ::Float64,η::Float64) = (1.,ξ,η,ξ^2,ξ*η,η^2,ξ^3,ξ^2*η,ξ*η^2,η^3)
@inline get∂𝒑₁∂ξ(ap::ReproducingKernel{:Quartic2D},ξ::SNode) = get∂𝒑₁∂ξ(ap,ξ.ξ,ξ.η)
@inline get∂𝒑₁∂ξ(::ReproducingKernel{:Quartic2D},ξ::Float64,η::Float64) = (0.,1.,0.,2.0*ξ,η,0.,3.0*ξ^2,2.0*ξ*η,η^2,0.)
@inline get∂𝒑₁∂η(ap::ReproducingKernel{:Quartic2D},ξ::SNode) = get∂𝒑₁∂η(ap,ξ.ξ,ξ.η)
@inline get∂𝒑₁∂η(::ReproducingKernel{:Quartic2D},ξ::Float64,η::Float64) = (0.,0.,1.,0.,ξ,2.0*η,0.,ξ^2,2.0*ξ*η,3.0*η^2)

@inline get𝑛𝒑₂(::ReproducingKernel{:Quartic2D}) = 6
@inline get𝒑₂(ap::ReproducingKernel{:Quartic2D},ξ::SNode) = get𝒑₂(ap,ξ.ξ,ξ.η)
@inline get𝒑₂(ap::ReproducingKernel{:Quartic2D},ξ::NTuple{3,Float64}) = get𝒑₂(ap,ξ[1],ξ[2])
@inline get𝒑₂(::ReproducingKernel{:Quartic2D},ξ::Float64,η::Float64) = (1.,ξ,η,ξ^2,ξ*η,η^2)
@inline get∂𝒑₂∂ξ(ap::ReproducingKernel{:Quartic2D},ξ::SNode) = get∂𝒑₂∂ξ(ap,ξ.ξ,ξ.η)
@inline get∂𝒑₂∂ξ(ap::ReproducingKernel{:Quartic2D},ξ::Float64,η::Float64) = (0.,1.,0.,2.0*ξ,η,0.)
@inline get∂𝒑₂∂η(ap::ReproducingKernel{:Quartic2D},ξ::SNode) = get∂𝒑₂∂η(ap,ξ.ξ,ξ.η)
@inline get∂𝒑₂∂η(ap::ReproducingKernel{:Quartic2D},ξ::Float64,η::Float64) = (0.,0.,1.,0.,ξ,2.0*η)
@inline get∂²𝒑₂∂ξ²(ap::ReproducingKernel{:Quartic2D},ξ::Any) = (0.,0.,0.,2.,0.,0.)
@inline get∂²𝒑₂∂ξ∂η(ap::ReproducingKernel{:Quartic2D},ξ::Any) = (0.,0.,0.,0.,1.,0.)
@inline get∂²𝒑₂∂η²(ap::ReproducingKernel{:Quartic2D},ξ::Any) = (0.,0.,0.,0.,0.,2.)

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

function get∇₂𝜙(ap::ReproducingKernel{𝒑,:□,𝜙},x::Node,Δx::NTuple{3,Float64}) where {𝒑,𝜙}
    rx = abs(Δx[1])/x.s₁
    ry = abs(Δx[2])/x.s₂
    ∂rx = sign(Δx[1])/x.s₁
    ∂ry = sign(Δx[2])/x.s₂
    wx = get𝜙ᵣ(ap,rx)
    wy = get𝜙ᵣ(ap,ry)
    ∂wx = get∂𝜙∂r(ap,rx)*∂rx
    ∂wy = get∂𝜙∂r(ap,ry)*∂ry
    return wx*wy, ∂wx*wy, wx*∂wy
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

function get∇²₂𝜙(ap::ReproducingKernel{𝒑,:□,𝜙},x::Node,Δx::NTuple{3,Float64}) where {𝒑,𝜙}
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
    return wx*wy, ∂wx*wy, wx*∂wy, ∂²wx*wy, ∂wx*∂wy, wx*∂²wy
end

function get∇³₂𝜙(ap::ReproducingKernel{𝒑,:□,𝜙},x::Node,Δx::NTuple{3,Float64}) where {𝒑,𝜙}
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

function get∇³𝜙(ap::ReproducingKernel{𝒑,:□,𝜙},x::Node,Δx::NTuple{3,Float64}) where {𝒑,𝜙}
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
    ∂²wx = get∂𝜙∂r(ap,rx)*∂rx^2
    ∂²wy = get∂𝜙∂r(ap,ry)*∂ry^2
    ∂³wx = get∂𝜙∂r(ap,rx)*∂rx^3
    ∂³wy = get∂𝜙∂r(ap,ry)*∂ry^3
    return wx*wy*wz, ∂wx*wy*wz, wx*∂wy*wz, ∂²wx*wy*wz, ∂wx*∂wy*wz, wx*∂²wy*wz, ∂³wx*wy*wz, ∂²wx*∂wy*wz, ∂wx*∂²wy*wz, wx*∂³wy*wz
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

function cal𝗠!(ap::ReproducingKernel,x::SNode)
    𝓒 = ap.𝓒
    𝗠 = get𝗠(ap,:𝗠)
    n = get𝑛𝒑(ap)
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
    return 𝗠
end

# function cal∇₁𝗠!(ap::ReproducingKernel,x::NTuple{3,Float64})
#     𝓒 = ap.𝓒
#     𝗠 = ap.𝗠[:∂1]
#     ∂𝗠∂x = ap.𝗠[:∂x]
#     n = get𝑛𝒑(ap)
#     fill!(𝗠,0.)
#     fill!(∂𝗠∂x,0.)
#     for xᵢ in 𝓒
#         Δx = x - xᵢ
#         𝒑, ∂𝒑∂x = get∇₁𝒑(ap,Δx)
#         𝜙, ∂𝜙∂x = get∇₁𝜙(ap,xᵢ,Δx)
#         for I in 1:n
#             for J in 1:I
#                 𝗠[I,J] += 𝜙*𝒑[I]*𝒑[J]
#                 ∂𝗠∂x[I,J] += ∂𝜙∂x*𝒑[I]*𝒑[J] + 𝜙*∂𝒑∂x[I]*𝒑[J] + 𝜙*𝒑[I]*∂𝒑∂x[J]
#             end
#         end
#     end
#     cholesky!(𝗠)
#     U = inverse!(𝗠)
#     ∂𝗠⁻¹∂x = - UUᵀAUUᵀ!(∂𝗠∂x,U)
#     𝗠⁻¹ = UUᵀ!(U)
#     return 𝗠⁻¹, ∂𝗠⁻¹∂x
# end

function cal∇₂𝗠!(ap::ReproducingKernel,x::SNode)
    𝓒 = ap.𝓒
    𝗠 = get𝗠(ap,:𝗠)
    ∂𝗠∂x = get𝗠(ap,:∂𝗠∂x)
    ∂𝗠∂y = get𝗠(ap,:∂𝗠∂y)
    n = get𝑛𝒑(ap)
    for xᵢ in 𝓒
        Δx = x - xᵢ
        𝒑, ∂𝒑∂x, ∂𝒑∂y = get∇𝒑(ap,Δx)
        𝜙, ∂𝜙∂x, ∂𝜙∂y = get∇𝜙(ap,xᵢ,Δx)
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

function cal∇𝗠!(ap::ReproducingKernel,x::SNode)
    𝓒 = ap.𝓒
    𝗠 = get𝗠(ap,:𝗠)
    ∂𝗠∂x = get𝗠(ap,:∂𝗠∂x)
    ∂𝗠∂y = get𝗠(ap,:∂𝗠∂y)
    ∂𝗠∂z = get𝗠(ap,:∂𝗠∂z)
    n = get𝑛𝒑(ap)
    for xᵢ in 𝓒
        Δx = x - xᵢ
        𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂𝒑∂z = get∇𝒑(ap,Δx)
        𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂𝜙∂z = get∇𝜙(ap,xᵢ,Δx)
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

function cal∇²₂𝗠!(ap::ReproducingKernel,x::SNode)
    𝓒 = ap.𝓒
    𝗠 = get𝗠(ap,:𝗠)
    ∂𝗠∂x = get𝗠(ap,:∂𝗠∂x)
    ∂𝗠∂y = get𝗠(ap,:∂𝗠∂y)
    ∂²𝗠∂x² = get𝗠(ap,:∂²𝝭∂x²)
    ∂²𝗠∂y² = get𝗠(ap,:∂²𝝭∂y²)
    ∂²𝗠∂x∂y = get𝗠(ap,:∂²𝝭∂x∂y)
    n = get𝑛𝒑(ap)
    for xᵢ in 𝓒
        Δx = x - xᵢ
        𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂²𝒑∂x², ∂²𝒑∂x∂y, ∂²𝒑∂y² = get∇²𝒑(ap,Δx)
        𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂²𝜙∂x², ∂²𝜙∂x∂y, ∂²𝜙∂y² = get∇²𝜙(ap,xᵢ,Δx)
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

function cal∇²𝗠!(ap::ReproducingKernel,x::SNode)
    𝓒 = ap.𝓒
    𝗠 = get𝗠(ap,:𝗠)
    ∂𝗠∂x = get𝗠(ap,:∂𝗠∂x)
    ∂𝗠∂y = get𝗠(ap,:∂𝗠∂y)
    ∂𝗠∂z = get𝗠(ap,:∂𝗠∂z)
    ∂²𝗠∂x² = get𝗠(ap,:∂²𝝭∂x²)
    ∂²𝗠∂y² = get𝗠(ap,:∂²𝝭∂y²)
    ∂²𝗠∂z² = get𝗠(ap,:∂²𝝭∂z²)
    ∂²𝗠∂x∂y = get𝗠(ap,:∂²𝝭∂x∂y)
    ∂²𝗠∂x∂z = get𝗠(ap,:∂²𝝭∂x∂z)
    ∂²𝗠∂y∂z = get𝗠(ap,:∂²𝝭∂y∂z)
    n = get𝑛𝒑(ap)
    for xᵢ in 𝓒
        Δx = x - xᵢ
        𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂²𝒑∂x², ∂²𝒑∂x∂y, ∂²𝒑∂y², ∂𝒑∂z, ∂²𝒑∂x∂z, ∂²𝒑∂y∂z, ∂²𝒑∂z² = get∇²𝒑(ap,Δx)
        𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂²𝜙∂x², ∂²𝜙∂x∂y, ∂²𝜙∂y², ∂𝜙∂z, ∂²𝜙∂x∂z, ∂²𝜙∂y∂z, ∂²𝜙∂z² = get∇²𝜙(ap,xᵢ,Δx)
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

function cal∇³𝗠!(ap::ReproducingKernel,x::SNode)
    𝓒 = ap.𝓒
    𝗠 = get𝗠(ap,:𝗠)
    ∂𝗠∂x = get𝗠(ap,:∂𝗠∂x)
    ∂𝗠∂y = get𝗠(ap,:∂𝗠∂y)
    ∂²𝗠∂x² = get𝗠(ap,:∂²𝝭∂x²)
    ∂²𝗠∂x∂y = get𝗠(ap,:∂²𝝭∂x∂y)
    ∂²𝗠∂y² = get𝗠(ap,:∂²𝝭∂y²)
    ∂³𝗠∂x³ = get𝗠(ap,:∂³𝝭∂x³)
    ∂³𝗠∂x²∂y = get𝗠(ap,:∂³𝝭∂x²∂y)
    ∂³𝗠∂x∂y² = get𝗠(ap,:∂³𝝭∂x∂y²)
    ∂³𝗠∂y³ = get𝗠(ap,:∂³𝝭∂y³)
    n = get𝑛𝒑(ap)
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

"""
cal𝗚!
"""
function cal𝗚!(ap::ReproducingKernel{:Linear1D,𝑠,𝜙,:Seg2}) where {𝑠,𝜙}
    𝗚⁻¹ = get𝗚(ap,:∇̃)
    𝐿 = ap.𝓖[1].𝐿
    𝗚⁻¹[1] =  1.0/𝐿
    return 𝗚⁻¹
end

function cal𝗚!(ap::ReproducingKernel{:Quadratic1D,𝑠,𝜙,:Seg2}) where {𝑠,𝜙}
    𝗚⁻¹ = get𝗚(ap,:∇̃)
    𝐿 = ap.𝓖[1].𝐿
    𝗚⁻¹[1] =  4.0/𝐿
    𝗚⁻¹[2] = -6.0/𝐿
    𝗚⁻¹[3] = 12.0/𝐿
    return 𝗚⁻¹
end

function cal𝗚!(ap::ReproducingKernel{:Cubic1D,𝑠,𝜙,:Seg2}) where {𝑠,𝜙}
    𝗚⁻¹ = get𝗚(ap,:∇̃)
    𝐿 = ap.𝓖[1].𝐿
    𝗚⁻¹[1] =    9.0/𝐿
    𝗚⁻¹[2] =  -36.0/𝐿
    𝗚⁻¹[3] =  192.0/𝐿
    𝗚⁻¹[4] =   30.0/𝐿
    𝗚⁻¹[5] = -180.0/𝐿
    𝗚⁻¹[6] =  180.0/𝐿
    return 𝗚⁻¹
end

function cal𝗚!(ap::ReproducingKernel{:Linear2D,𝑠,𝜙,:Tri3}) where {𝑠,𝜙}
    𝗚⁻¹ = get𝗚(ap,:∇̃)
    𝐴 = ap.𝓖[1].𝐴
    𝗚⁻¹[1] = 1.0/𝐴
    return 𝗚⁻¹
end

function cal𝗚!(ap::ReproducingKernel{:Quadratic2D,𝑠,𝜙,:Tri3}) where {𝑠,𝜙}
    𝗚⁻¹ = get𝗚(ap,:∇̃)
    𝐴 = ap.𝓖[1].𝐴
    𝗚⁻¹[1] =   9.0/𝐴
    𝗚⁻¹[2] = -12.0/𝐴
    𝗚⁻¹[3] =  24.0/𝐴
    𝗚⁻¹[4] = -12.0/𝐴
    𝗚⁻¹[5] =  12.0/𝐴
    𝗚⁻¹[6] =  24.0/𝐴
    return 𝗚⁻¹
end

function cal𝗚!(ap::ReproducingKernel{:Cubic2D,𝑠,𝜙,:Tri3}) where {𝑠,𝜙}
    𝗚⁻¹ = get𝗚(ap,:∇̃)
    𝐴 = ap.𝓖[1].𝐴
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

function cal𝗚₂!(ap::ReproducingKernel{:Quadratic2D,𝑠,𝜙,:Tri3}) where {𝑠,𝜙}
    𝗚⁻¹ = get𝗚(ap,:∇̃²)
    𝐴 = ap.𝓖[1].𝐴
    𝗚⁻¹[1] = 1.0/𝐴
    return 𝗚⁻¹
end

function cal𝗚₂!(ap::ReproducingKernel{:Cubic2D,𝑠,𝜙,:Tri3}) where {𝑠,𝜙}
    𝗚⁻¹ = get𝗚(ap,:∇̃²)
    𝐴 = ap.𝓖[1].𝐴
    𝗚⁻¹[1] =   9.0/𝐴
    𝗚⁻¹[2] = -12.0/𝐴
    𝗚⁻¹[3] =  24.0/𝐴
    𝗚⁻¹[4] = -12.0/𝐴
    𝗚⁻¹[5] =  12.0/𝐴
    𝗚⁻¹[6] =  24.0/𝐴
    return 𝗚⁻¹
end

function cal∇𝗚₂!(ap::ReproducingKernel{:Cubic2D,𝑠,𝜙,:Tri3}) where {𝑠,𝜙}
    𝗚⁻¹ = get𝗚(ap,:∇̃²)
    𝗚⁻¹∂ξ = get𝗚(ap,:∂∇̃²∂ξ)
    𝗚⁻¹∂η = get𝗚(ap,:∂∇̃²∂η)
    𝐴 = ap.𝓖[1].𝐴
    𝗚⁻¹[1] =   9.0/𝐴
    𝗚⁻¹[2] = -12.0/𝐴
    𝗚⁻¹[3] =  24.0/𝐴
    𝗚⁻¹[4] = -12.0/𝐴
    𝗚⁻¹[5] =  12.0/𝐴
    𝗚⁻¹[6] =  24.0/𝐴

    𝗚⁻¹∂ξ[1] =   9.0/𝐴
    𝗚⁻¹∂ξ[2] = -12.0/𝐴
    𝗚⁻¹∂ξ[3] =  24.0/𝐴
    𝗚⁻¹∂ξ[4] = -12.0/𝐴
    𝗚⁻¹∂ξ[5] =  12.0/𝐴
    𝗚⁻¹∂ξ[6] =  24.0/𝐴

    𝗚⁻¹∂η[1] =   9.0/𝐴
    𝗚⁻¹∂η[2] = -12.0/𝐴
    𝗚⁻¹∂η[3] =  24.0/𝐴
    𝗚⁻¹∂η[4] = -12.0/𝐴
    𝗚⁻¹∂η[5] =  12.0/𝐴
    𝗚⁻¹∂η[6] =  24.0/𝐴

    return 𝗚⁻¹,𝗚⁻¹∂ξ,𝗚⁻¹∂η
end

function cal𝗚₂!(ap::ReproducingKernel{:Quartic2D,𝑠,𝜙,:Tri3}) where {𝑠,𝜙}
    𝗚⁻¹ = get𝗚(ap,:∇̃²)
    𝐴 = ap.𝓖[1].𝐴
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

function cal∇𝗚₂!(ap::ReproducingKernel{:Quartic2D,𝑠,𝜙,:Tri3}) where {𝑠,𝜙}
    𝗚⁻¹ = get𝗚(ap,:∇̃²)
    𝗚⁻¹∂ξ = get𝗚(ap,:∂∇̃²∂ξ)
    𝗚⁻¹∂η = get𝗚(ap,:∂∇̃²∂η)
    𝐴 = ap.𝓖[1].𝐴
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

    𝗚⁻¹∂ξ[1] =   36.0/𝐴
    𝗚⁻¹∂ξ[2] = -120.0/𝐴
    𝗚⁻¹∂ξ[3] =  600.0/𝐴
    𝗚⁻¹∂ξ[4] = -120.0/𝐴
    𝗚⁻¹∂ξ[5] =  300.0/𝐴
    𝗚⁻¹∂ξ[6] =  600.0/𝐴
    𝗚⁻¹∂ξ[7] =   90.0/𝐴
    𝗚⁻¹∂ξ[8] = -540.0/𝐴
    𝗚⁻¹∂ξ[9] = -180.0/𝐴
    𝗚⁻¹∂ξ[10] =  540.0/𝐴
    𝗚⁻¹∂ξ[11] =  180.0/𝐴
    𝗚⁻¹∂ξ[12] = -720.0/𝐴
    𝗚⁻¹∂ξ[13] = -720.0/𝐴
    𝗚⁻¹∂ξ[14] =  540.0/𝐴
    𝗚⁻¹∂ξ[15] = 1440.0/𝐴
    𝗚⁻¹∂ξ[16] =   90.0/𝐴
    𝗚⁻¹∂ξ[17] = -180.0/𝐴
    𝗚⁻¹∂ξ[18] = -540.0/𝐴
    𝗚⁻¹∂ξ[19] =   90.0/𝐴
    𝗚⁻¹∂ξ[20] =  540.0/𝐴
    𝗚⁻¹∂ξ[21] =  540.0/𝐴

    𝗚⁻¹∂η[1] =   36.0/𝐴
    𝗚⁻¹∂η[2] = -120.0/𝐴
    𝗚⁻¹∂η[3] =  600.0/𝐴
    𝗚⁻¹∂η[4] = -120.0/𝐴
    𝗚⁻¹∂η[5] =  300.0/𝐴
    𝗚⁻¹∂η[6] =  600.0/𝐴
    𝗚⁻¹∂η[7] =   90.0/𝐴
    𝗚⁻¹∂η[8] = -540.0/𝐴
    𝗚⁻¹∂η[9] = -180.0/𝐴
    𝗚⁻¹∂η[10] =  540.0/𝐴
    𝗚⁻¹∂η[11] =  180.0/𝐴
    𝗚⁻¹∂η[12] = -720.0/𝐴
    𝗚⁻¹∂η[13] = -720.0/𝐴
    𝗚⁻¹∂η[14] =  540.0/𝐴
    𝗚⁻¹∂η[15] = 1440.0/𝐴
    𝗚⁻¹∂η[16] =   90.0/𝐴
    𝗚⁻¹∂η[17] = -180.0/𝐴
    𝗚⁻¹∂η[18] = -540.0/𝐴
    𝗚⁻¹∂η[19] =   90.0/𝐴
    𝗚⁻¹∂η[20] =  540.0/𝐴
    𝗚⁻¹∂η[21] =  540.0/𝐴

    return 𝗚⁻¹,𝗚⁻¹∂ξ,𝗚⁻¹∂η
end

"""
set𝝭!
"""
function set𝝭!(ap::ReproducingKernel,𝒙::SNode)
    𝓒 = ap.𝓒
    𝝭 = 𝒙[:𝝭]
    𝒑₀ᵀ𝗠⁻¹ = cal𝗠!(ap,𝒙)
    for (i,𝒙ᵢ) in enumerate(𝓒)
        Δ𝒙 = 𝒙 - 𝒙ᵢ
        𝒑= get𝒑(ap,Δ𝒙)
        𝜙 = get𝜙(ap,𝒙ᵢ,Δ𝒙)
        𝝭[i] = 𝒑₀ᵀ𝗠⁻¹*𝒑*𝜙
    end
end

"""
set∇𝝭!
"""
function set∇𝝭!(ap::ReproducingKernel,𝒙::SNode)
    𝓒 = ap.𝓒
    𝝭 = 𝒙[:𝝭]
    ∂𝝭∂x = 𝒙[:∂𝝭∂x]
    ∂𝝭∂y = 𝒙[:∂𝝭∂y]
    ∂𝝭∂z = 𝒙[:∂𝝭∂z]
    𝒑₀ᵀ𝗠⁻¹, 𝒑₀ᵀ∂𝗠⁻¹∂x, 𝒑₀ᵀ∂𝗠⁻¹∂y, 𝒑₀ᵀ∂𝗠⁻¹∂z= cal∇𝗠!(ap,𝒙)
    for (i,𝒙ᵢ) in enumerate(𝓒)
        Δ𝒙 = 𝒙 - 𝒙ᵢ
        𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂𝒑∂z = get∇𝒑(ap,Δ𝒙)
        𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂𝜙∂z = get∇𝜙(ap,𝒙ᵢ,Δ𝒙)
        𝝭[i] = 𝒑₀ᵀ𝗠⁻¹*𝒑*𝜙
        ∂𝝭∂x[i] = 𝒑₀ᵀ∂𝗠⁻¹∂x*𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂x*𝜙 + 𝒑₀ᵀ𝗠⁻¹*𝒑*∂𝜙∂x
        ∂𝝭∂y[i] = 𝒑₀ᵀ∂𝗠⁻¹∂y*𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂y*𝜙 + 𝒑₀ᵀ𝗠⁻¹*𝒑*∂𝜙∂y
        ∂𝝭∂z[i] = 𝒑₀ᵀ∂𝗠⁻¹∂z*𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂z*𝜙 + 𝒑₀ᵀ𝗠⁻¹*𝒑*∂𝜙∂z
    end
end

function set∇₂𝝭!(ap::ReproducingKernel,𝒙::SNode)
    𝓒 = ap.𝓒
    𝝭 = 𝒙[:𝝭]
    ∂𝝭∂x = 𝒙[:∂𝝭∂x]
    ∂𝝭∂y = 𝒙[:∂𝝭∂y]
    𝒑₀ᵀ𝗠⁻¹, 𝒑₀ᵀ∂𝗠⁻¹∂x, 𝒑₀ᵀ∂𝗠⁻¹∂y = cal∇₂𝗠!(ap,𝒙)
    for (i,𝒙ᵢ) in enumerate(𝓒)
        Δ𝒙 = 𝒙 - 𝒙ᵢ
        𝒑, ∂𝒑∂x, ∂𝒑∂y = get∇𝒑(ap,Δ𝒙)
        𝜙, ∂𝜙∂x, ∂𝜙∂y = get∇𝜙(ap,𝒙ᵢ,Δ𝒙)
        𝝭[i] = 𝒑₀ᵀ𝗠⁻¹*𝒑*𝜙
        ∂𝝭∂x[i] = 𝒑₀ᵀ∂𝗠⁻¹∂x*𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂x*𝜙 + 𝒑₀ᵀ𝗠⁻¹*𝒑*∂𝜙∂x
        ∂𝝭∂y[i] = 𝒑₀ᵀ∂𝗠⁻¹∂y*𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂y*𝜙 + 𝒑₀ᵀ𝗠⁻¹*𝒑*∂𝜙∂y
    end
end

"""
set∇²𝝭!
"""
function set∇²𝝭!(ap::ReproducingKernel,𝒙::SNode)
    𝓒 = ap.𝓒
    𝝭 = 𝒙[:𝝭]
    ∂𝝭∂x = 𝒙[:∂𝝭∂x]
    ∂𝝭∂y = 𝒙[:∂𝝭∂y]
    ∂𝝭∂z = 𝒙[:∂𝝭∂z]
    ∂²𝝭∂x² = 𝒙[:∂²𝝭∂x²]
    ∂²𝝭∂y² = 𝒙[:∂²𝝭∂y²]
    ∂²𝝭∂z² = 𝒙[:∂²𝝭∂z²]
    ∂²𝝭∂x∂y = 𝒙[:∂²𝝭∂x∂y]
    ∂²𝝭∂x∂z = 𝒙[:∂²𝝭∂x∂z]
    ∂²𝝭∂y∂z = 𝒙[:∂²𝝭∂y∂z]
    𝒑₀ᵀ𝗠⁻¹, 𝒑₀ᵀ∂𝗠⁻¹∂x, 𝒑₀ᵀ∂𝗠⁻¹∂y, 𝒑₀ᵀ∂²𝗠⁻¹∂x², 𝒑₀ᵀ∂²𝗠⁻¹∂x∂y, 𝒑₀ᵀ∂²𝗠⁻¹∂y², 𝒑₀ᵀ∂𝗠⁻¹∂z, 𝒑₀ᵀ∂²𝗠⁻¹∂x∂z, 𝒑₀ᵀ∂²𝗠⁻¹∂y∂z, 𝒑₀ᵀ∂²𝗠⁻¹∂z² = cal∇²𝗠!(ap,𝒙)
    for (i,𝒙ᵢ) in enumerate(𝓒)
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

        ∂²𝝭∂x∂z[i] = 𝒑₀ᵀ∂²𝗠⁻¹∂x∂z𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x∂z*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂²𝜙∂x∂z + 𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂z*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂z∂𝒑∂x*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*∂𝜙∂z + 𝒑₀ᵀ∂𝗠⁻¹∂z𝒑*∂𝜙∂x + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*∂𝜙∂z +𝒑₀ᵀ𝗠⁻¹∂𝒑∂z*∂𝜙∂x

        ∂²𝝭∂y∂z[i] = 𝒑₀ᵀ∂²𝗠⁻¹∂y∂z𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂y∂z*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂²𝜙∂y∂z + 𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂z*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂z∂𝒑∂y*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*∂𝜙∂z + 𝒑₀ᵀ∂𝗠⁻¹∂z𝒑*∂𝜙∂y + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*∂𝜙∂z +𝒑₀ᵀ𝗠⁻¹∂𝒑∂z*∂𝜙∂y
    end
end

function set∇²₂𝝭!(ap::ReproducingKernel,𝒙::SNode)
    𝓒 = ap.𝓒
    𝝭 = 𝒙[:𝝭]
    ∂𝝭∂x = 𝒙[:∂𝝭∂x]
    ∂𝝭∂y = 𝒙[:∂𝝭∂y]
    ∂²𝝭∂x² = 𝒙[:∂²𝝭∂x²]
    ∂²𝝭∂y² = 𝒙[:∂²𝝭∂y²]
    ∂²𝝭∂x∂y = 𝒙[:∂²𝝭∂x∂y]
    𝒑₀ᵀ𝗠⁻¹, 𝒑₀ᵀ∂𝗠⁻¹∂x, 𝒑₀ᵀ∂𝗠⁻¹∂y, 𝒑₀ᵀ∂²𝗠⁻¹∂x², 𝒑₀ᵀ∂²𝗠⁻¹∂x∂y, 𝒑₀ᵀ∂²𝗠⁻¹∂y² = cal∇²₂𝗠!(ap,𝒙)
    for (i,𝒙ᵢ) in enumerate(𝓒)
        Δ𝒙 = 𝒙 - 𝒙ᵢ
        𝒑, ∂𝒑∂x, ∂𝒑∂y, ∂²𝒑∂x², ∂²𝒑∂x∂y, ∂²𝒑∂y² = get∇²𝒑(ap,Δ𝒙)
        𝜙, ∂𝜙∂x, ∂𝜙∂y, ∂²𝜙∂x², ∂²𝜙∂x∂y, ∂²𝜙∂y² = get∇²𝜙(ap,𝒙ᵢ,Δ𝒙)
        𝒑₀ᵀ𝗠⁻¹𝒑 = 𝒑₀ᵀ𝗠⁻¹*𝒑
        𝒑₀ᵀ∂𝗠⁻¹∂x𝒑 = 𝒑₀ᵀ∂𝗠⁻¹∂x*𝒑
        𝒑₀ᵀ∂𝗠⁻¹∂y𝒑 = 𝒑₀ᵀ∂𝗠⁻¹∂y*𝒑
        𝒑₀ᵀ𝗠⁻¹∂𝒑∂x = 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂x
        𝒑₀ᵀ𝗠⁻¹∂𝒑∂y = 𝒑₀ᵀ𝗠⁻¹*∂𝒑∂y
        𝒑₀ᵀ∂²𝗠⁻¹∂x²𝒑 = 𝒑₀ᵀ∂²𝗠⁻¹∂x²*𝒑
        𝒑₀ᵀ∂²𝗠⁻¹∂y²𝒑 = 𝒑₀ᵀ∂²𝗠⁻¹∂y²*𝒑
        𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x² = 𝒑₀ᵀ𝗠⁻¹*∂²𝒑∂x²
        𝒑₀ᵀ𝗠⁻¹∂²𝒑∂y² = 𝒑₀ᵀ𝗠⁻¹*∂²𝒑∂y²
        𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂x = 𝒑₀ᵀ∂𝗠⁻¹∂x*∂𝒑∂x
        𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂y = 𝒑₀ᵀ∂𝗠⁻¹∂y*∂𝒑∂y
        𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂y = 𝒑₀ᵀ∂𝗠⁻¹∂x*∂𝒑∂y
        𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂x = 𝒑₀ᵀ∂𝗠⁻¹∂y*∂𝒑∂x
        𝒑₀ᵀ∂²𝗠⁻¹∂x∂y𝒑 = 𝒑₀ᵀ∂²𝗠⁻¹∂x∂y*𝒑
        𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x∂y = 𝒑₀ᵀ𝗠⁻¹*∂²𝒑∂x∂y

        𝝭[i] = 𝒑₀ᵀ𝗠⁻¹𝒑*𝜙
        ∂𝝭∂x[i] = 𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂𝜙∂x
        ∂𝝭∂y[i] = 𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂𝜙∂y

        ∂²𝝭∂x²[i] = 𝒑₀ᵀ∂²𝗠⁻¹∂x²𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x²*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂²𝜙∂x² + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂x*𝜙 + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*∂𝜙∂x + 2.0*𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*∂𝜙∂x

        ∂²𝝭∂y²[i] = 𝒑₀ᵀ∂²𝗠⁻¹∂y²𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂y²*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂²𝜙∂y² + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂y*𝜙 + 2.0*𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*∂𝜙∂y + 2.0*𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*∂𝜙∂y

        ∂²𝝭∂x∂y[i] = 𝒑₀ᵀ∂²𝗠⁻¹∂x∂y𝒑*𝜙 + 𝒑₀ᵀ𝗠⁻¹∂²𝒑∂x∂y*𝜙 + 𝒑₀ᵀ𝗠⁻¹𝒑*∂²𝜙∂x∂y + 𝒑₀ᵀ∂𝗠⁻¹∂x∂𝒑∂y*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂y∂𝒑∂x*𝜙 + 𝒑₀ᵀ∂𝗠⁻¹∂x𝒑*∂𝜙∂y + 𝒑₀ᵀ∂𝗠⁻¹∂y𝒑*∂𝜙∂x + 𝒑₀ᵀ𝗠⁻¹∂𝒑∂x*∂𝜙∂y +𝒑₀ᵀ𝗠⁻¹∂𝒑∂y*∂𝜙∂x
    end
end

function set∇³𝝭!(ap::ReproducingKernel,𝒙::SNode)
    𝓒 = ap.𝓒
    𝝭 = 𝒙[:𝝭]
    ∂𝝭∂x = 𝒙[:∂𝝭∂x]
    ∂𝝭∂y = 𝒙[:∂𝝭∂y]
    ∂²𝝭∂x² = 𝒙[:∂²𝝭∂x²]
    ∂²𝝭∂y² = 𝒙[:∂²𝝭∂y²]
    ∂²𝝭∂x∂y = 𝒙[:∂²𝝭∂x∂y]
    ∂³𝝭∂x³ = 𝒙[:∂³𝝭∂x³]
    ∂³𝝭∂x²∂y = 𝒙[:∂³𝝭∂x²∂y]
    ∂³𝝭∂x∂y² = 𝒙[:∂³𝝭∂x∂y²]
    ∂³𝝭∂y³ = 𝒙[:∂³𝝭∂y³]
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
end

function set∇̂³𝝭!(ap::ReproducingKernel,𝒙::SNode)
    𝓒 = ap.𝓒
    𝝭 = 𝒙[:𝝭]
    ∂𝝭∂x = 𝒙[:∂𝝭∂x]
    ∂𝝭∂y = 𝒙[:∂𝝭∂y]
    ∂²𝝭∂x² = 𝒙[:∂²𝝭∂x²]
    ∂²𝝭∂y² = 𝒙[:∂²𝝭∂y²]
    ∂²𝝭∂x∂y = 𝒙[:∂²𝝭∂x∂y]
    ∂³𝝭∂x³ = 𝒙[:∂³𝝭∂x³]
    ∂³𝝭∂x²∂y = 𝒙[:∂³𝝭∂x²∂y]
    ∂³𝝭∂x∂y² = 𝒙[:∂³𝝭∂x∂y²]
    ∂³𝝭∂y³ = 𝒙[:∂³𝝭∂y³]

    n = get𝑛𝒑(ap)
    𝒑₀ᵀ𝗠⁻¹ = cal𝗠!(ap,𝒙)
    for (i,𝒙ᵢ) in enumerate(𝓒)
        Δ𝒙 = 𝒙 - 𝒙ᵢ
        𝒑= get𝒑(ap,Δ𝒙)
        𝜙 = get𝜙(ap,𝒙ᵢ,Δ𝒙)
        ∂𝝭∂x_ = 0.0
        ∂𝝭∂y_ = 0.0
        ∂²𝝭∂x²_ = 0.0
        ∂²𝝭∂x∂y_ = 0.0
        ∂²𝝭∂y²_ = 0.0
        ∂³𝝭∂x³_ = 0.0
        ∂³𝝭∂x²∂y_ = 0.0
        ∂³𝝭∂x∂y²_ = 0.0
        ∂³𝝭∂y³_ = 0.0
        for j in 1:n
            ∂𝝭∂x_ -= 𝒑₀ᵀ𝗠⁻¹[2,j]*𝒑[j]*𝜙
            ∂𝝭∂y_ -= 𝒑₀ᵀ𝗠⁻¹[3,j]*𝒑[j]*𝜙
            ∂²𝝭∂x²_ += 2*𝒑₀ᵀ𝗠⁻¹[4,j]*𝒑[j]*𝜙
            ∂²𝝭∂x∂y_ += 𝒑₀ᵀ𝗠⁻¹[5,j]*𝒑[j]*𝜙
            ∂²𝝭∂y²_ += 2*𝒑₀ᵀ𝗠⁻¹[6,j]*𝒑[j]*𝜙
            ∂³𝝭∂x³_ -= 6*𝒑₀ᵀ𝗠⁻¹[7,j]*𝒑[j]*𝜙
            ∂³𝝭∂x²∂y_ -= 2*𝒑₀ᵀ𝗠⁻¹[8,j]*𝒑[j]*𝜙
            ∂³𝝭∂x∂y²_ -= 2*𝒑₀ᵀ𝗠⁻¹[9,j]*𝒑[j]*𝜙
            ∂³𝝭∂y³_ -= 6*𝒑₀ᵀ𝗠⁻¹[10,j]*𝒑[j]*𝜙
        end
        𝝭[i] = 𝒑₀ᵀ𝗠⁻¹*𝒑*𝜙
        ∂𝝭∂x[i] = ∂𝝭∂x_
        ∂𝝭∂y[i] = ∂𝝭∂y_
        ∂²𝝭∂x²[i] = ∂²𝝭∂x²_
        ∂²𝝭∂x∂y[i] = ∂²𝝭∂x∂y_
        ∂²𝝭∂y²[i] = ∂²𝝭∂y²_
        ∂³𝝭∂x³[i] = ∂³𝝭∂x³_
        ∂³𝝭∂x²∂y[i] = ∂³𝝭∂x²∂y_
        ∂³𝝭∂x∂y²[i] = ∂³𝝭∂x∂y²_
        ∂³𝝭∂y³[i] = ∂³𝝭∂y³_
    end
end


"""
set∇̃𝝭!
"""
function set∇̃𝝭!(gp::ReproducingKernel{𝒑,𝑠,𝜙,:Seg2},ap::ReproducingKernel{𝒑,𝑠,𝜙,:Seg2}) where {𝒑,𝑠,𝜙}
    𝓒 = gp.𝓒
    𝓖 = gp.𝓖
    for ξ̂ in 𝓖
        𝒒̂ = get𝒑₁(gp,ξ̂)
        𝗚⁻¹ = cal𝗚!(gp)
        𝒒̂ᵀ𝗚⁻¹ = 𝒒̂*𝗚⁻¹
        ∂𝝭∂x = ξ̂[:∂𝝭∂x]
        for i in 1:length(𝓒)
            ∂𝝭∂x[i] = 0.0
        end
        for ξ in ap.𝓖
            w = ξ.w/2
            wᵇ = ξ.wᵇ
            D₁ = ξ.D₁
            𝝭 = ξ[:𝝭]
            𝒒, ∂𝒒∂ξ = get∇𝒑₁(gp,ξ)
            W₁ = 𝒒̂ᵀ𝗚⁻¹*𝒒*D₁*wᵇ + 𝒒̂ᵀ𝗚⁻¹*∂𝒒∂ξ*n₁*w
            for i in 1:length(𝓒)
                ∂𝝭∂x[i] += 𝝭[i]*W₁
            end
        end
    end
end

function set∇̃𝝭!(gp::ReproducingKernel{𝒑,𝑠,𝜙,:Tri3},ap::ReproducingKernel{𝒑,𝑠,𝜙,:Tri3}) where {𝒑,𝑠,𝜙}
    𝓒 = gp.𝓒
    𝓖 = gp.𝓖
    for ξ̂ in 𝓖
        𝒒̂ = get𝒑₁(gp,ξ̂)
        𝗚⁻¹ = cal𝗚!(gp)
        𝒒̂ᵀ𝗚⁻¹ = 𝒒̂*𝗚⁻¹
        ∂𝝭∂x = ξ̂[:∂𝝭∂x]
        ∂𝝭∂y = ξ̂[:∂𝝭∂y]
        for i in 1:length(𝓒)
            ∂𝝭∂x[i] = 0.0
            ∂𝝭∂y[i] = 0.0
        end
        for ξ in ap.𝓖
            w = ξ.w
            wᵇ = ξ.wᵇ
            𝝭 = ξ[:𝝭]
            𝒒, ∂𝒒∂ξ, ∂𝒒∂η = get∇𝒑₁(ap,ξ)
            𝒒̂ᵀ𝗚⁻¹𝒒 =  𝒒̂ᵀ𝗚⁻¹*𝒒
            𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ = 𝒒̂ᵀ𝗚⁻¹*∂𝒒∂ξ
            𝒒̂ᵀ𝗚⁻¹∂𝒒∂η = 𝒒̂ᵀ𝗚⁻¹*∂𝒒∂η
            D₁ = ξ.D₁
            D₂ = ξ.D₂
            D₁₁ = ξ.D₁₁
            D₂₁ = ξ.D₂₁
            D₁₂ = ξ.D₁₂
            D₂₂ = ξ.D₂₂
            b₁ = 𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ*D₁₁ + 𝒒̂ᵀ𝗚⁻¹∂𝒒∂η*D₂₁
            b₂ = 𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ*D₁₂ + 𝒒̂ᵀ𝗚⁻¹∂𝒒∂η*D₂₂
            W₁ = 𝒒̂ᵀ𝗚⁻¹𝒒*D₁*wᵇ + b₁*w/2
            W₂ = 𝒒̂ᵀ𝗚⁻¹𝒒*D₂*wᵇ + b₂*w/2
            for i in 1:length(𝓒)
                ∂𝝭∂x[i] += 𝝭[i]*W₁
                ∂𝝭∂y[i] += 𝝭[i]*W₂
            end
        end
    end
end

function set∇̃²𝝭!(gp::ReproducingKernel{𝒑,𝑠,𝜙,:Tri3},ap::ReproducingKernel{𝒑,𝑠,𝜙,:Tri3}) where {𝒑,𝑠,𝜙}
    x₁ = ap.𝓒[1].x;y₁ = ap.𝓒[1].y
    x₂ = ap.𝓒[2].x;y₂ = ap.𝓒[2].y
    x₃ = ap.𝓒[3].x;y₃ = ap.𝓒[3].y
    𝐴 = get𝐴(ap)
    n₁₁ = y₃-y₂;n₂₁ = y₁-y₃;n₃₁ = y₂-y₁
    n₁₂ = x₂-x₃;n₂₂ = x₃-x₁;n₃₂ = x₁-x₂
    s₁₁ = -n₁₂;s₂₁ = -n₂₂;s₃₁ = -n₃₂
    s₁₂ =  n₁₁;s₂₂ =  n₂₁;s₃₂ =  n₃₁
    𝐿₁² = n₁₁^2+n₁₂^2
    𝐿₂² = n₂₁^2+n₂₂^2
    𝐿₃² = n₃₁^2+n₃₂^2
    𝓒 = gp.𝓒
    𝓖 = gp.𝓖
    for ξ̂ in 𝓖
        𝒒̂ = get𝒑₂(gp,ξ̂)
        𝗚⁻¹ = cal𝗚₂!(gp)
        𝒒̂ᵀ𝗚⁻¹ = 𝒒̂*𝗚⁻¹

        ∂²𝝭∂x² = ξ̂[:∂²𝝭∂x²]
        ∂²𝝭∂x∂y = ξ̂[:∂²𝝭∂x∂y]
        ∂²𝝭∂y² = ξ̂[:∂²𝝭∂y²]
        for i in 1:length(𝓒)
            ∂²𝝭∂x²[i] = 0.0
            ∂²𝝭∂x∂y[i] = 0.0
            ∂²𝝭∂y²[i] = 0.0
        end
        for ξ in ap.𝓖
            w = ξ.w
            wᵇ = ξ.wᵇ
            𝝭 = ξ[:𝝭]
            ∂𝝭∂x = ξ[:∂𝝭∂x]
            ∂𝝭∂y = ξ[:∂𝝭∂y]
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
            end
            if ξ.ξ == 1.0
                Δms₁₁ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₂₁*s₂₁/𝐿₂²-n₃₁*s₃₁/𝐿₃²)
                Δms₂₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₂₂*s₂₂/𝐿₂²-n₃₂*s₃₂/𝐿₃²)
                Δms₁₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*0.5*((n₂₁*s₂₂+n₂₂*s₂₁)/𝐿₂²-(n₃₁*s₃₂+n₃₂*s₃₁)/𝐿₃²)
            end
            if ξ.η == 1.0
                Δms₁₁ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₃₁*s₃₁/𝐿₃²-n₁₁*s₁₁/𝐿₁²)
                Δms₂₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₃₂*s₃₂/𝐿₃²-n₁₂*s₁₂/𝐿₁²)
                Δms₁₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*0.5*((n₃₁*s₃₂+n₃₂*s₃₁)/𝐿₃²-(n₁₁*s₁₂+n₁₂*s₁₁)/𝐿₁²)
            end
            if ξ.ξ+ξ.η ≈ 0.0
                Δms₁₁ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₁₁*s₁₁/𝐿₁²-n₂₁*s₂₁/𝐿₂²)
                Δms₂₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₁₂*s₁₂/𝐿₁²-n₂₂*s₂₂/𝐿₂²)
                Δms₁₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*0.5*((n₁₁*s₁₂+n₁₂*s₁₁)/𝐿₁²-(n₂₁*s₂₂+n₂₂*s₂₁)/𝐿₂²)
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
    end
end


function set∇∇̃²𝝭!(gp::ReproducingKernel{𝒑,𝑠,𝜙,:Tri3},ap::ReproducingKernel{𝒑,𝑠,𝜙,:Tri3}) where {𝒑,𝑠,𝜙}
    x₁ = ap.𝓒[1].x;y₁ = ap.𝓒[1].y
    x₂ = ap.𝓒[2].x;y₂ = ap.𝓒[2].y
    x₃ = ap.𝓒[3].x;y₃ = ap.𝓒[3].y
    𝐴 = get𝐴(ap)
    n₁₁ = y₃-y₂;n₂₁ = y₁-y₃;n₃₁ = y₂-y₁
    n₁₂ = x₂-x₃;n₂₂ = x₃-x₁;n₃₂ = x₁-x₂
    s₁₁ = -n₁₂;s₂₁ = -n₂₂;s₃₁ = -n₃₂
    s₁₂ =  n₁₁;s₂₂ =  n₂₁;s₃₂ =  n₃₁
    𝐿₁² = n₁₁^2+n₁₂^2
    𝐿₂² = n₂₁^2+n₂₂^2
    𝐿₃² = n₃₁^2+n₃₂^2
    𝓒 = gp.𝓒
    𝓖 = gp.𝓖
    for ξ̂ in 𝓖
        𝒒̂,∂𝒒̂∂ξ,∂𝒒̂∂η = get∇𝒑₂(gp,ξ̂)
        𝗚⁻¹,𝗚⁻¹∂ξ,𝗚⁻¹∂η = cal∇𝗚₂!(gp)
        𝒒̂ᵀ𝗚⁻¹ = 𝒒̂*𝗚⁻¹
        ∂𝒒̂ᵀ∂ξ𝗚⁻¹ = ∂𝒒̂∂ξ*𝗚⁻¹∂ξ
        ∂𝒒̂ᵀ∂η𝗚⁻¹ = ∂𝒒̂∂η*𝗚⁻¹∂η

        ∂²𝝭∂x² = ξ̂[:∂²𝝭∂x²]
        ∂²𝝭∂x∂y = ξ̂[:∂²𝝭∂x∂y]
        ∂²𝝭∂y² = ξ̂[:∂²𝝭∂y²]
        ∂∂²𝝭∂x²∂x = ξ̂[:∂∂²𝝭∂x²∂x]
        ∂∂²𝝭∂x²∂y = ξ̂[:∂∂²𝝭∂x²∂y]
        ∂∂²𝝭∂x∂y∂x = ξ̂[:∂∂²𝝭∂x∂y∂x]
        ∂∂²𝝭∂x∂y∂y = ξ̂[:∂∂²𝝭∂x∂y∂y]
        ∂∂²𝝭∂y²∂x = ξ̂[:∂∂²𝝭∂y²∂x]
        ∂∂²𝝭∂y²∂y = ξ̂[:∂∂²𝝭∂y²∂y]
        for i in 1:length(𝓒)
            ∂²𝝭∂x²[i] = 0.0
            ∂²𝝭∂x∂y[i] = 0.0
            ∂²𝝭∂y²[i] = 0.0
            ∂∂²𝝭∂x²∂x[i] = 0.0
            ∂∂²𝝭∂x²∂y[i] = 0.0
            ∂∂²𝝭∂x∂y∂x[i] = 0.0
            ∂∂²𝝭∂x∂y∂y[i] = 0.0
            ∂∂²𝝭∂y²∂x[i] = 0.0
            ∂∂²𝝭∂y²∂y[i] = 0.0
        end

        for ξ in ap.𝓖
            w = ξ.w
            wᵇ = ξ.wᵇ
            𝝭 = ξ[:𝝭]
            ∂𝝭∂x = ξ[:∂𝝭∂x]
            ∂𝝭∂y = ξ[:∂𝝭∂y]
            𝒒, ∂𝒒∂ξ, ∂𝒒∂η, ∂²𝒒∂ξ², ∂²𝒒∂ξ∂η, ∂²𝒒∂η² = get∇²𝒑₂(ap,ξ)

            𝒒̂ᵀ𝗚⁻¹𝒒 =  𝒒̂ᵀ𝗚⁻¹*𝒒
            𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ = 𝒒̂ᵀ𝗚⁻¹*∂𝒒∂ξ
            𝒒̂ᵀ𝗚⁻¹∂𝒒∂η = 𝒒̂ᵀ𝗚⁻¹*∂𝒒∂η
            𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ² = 𝒒̂ᵀ𝗚⁻¹*∂²𝒒∂ξ²
            𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ∂η = 𝒒̂ᵀ𝗚⁻¹*∂²𝒒∂ξ∂η
            𝒒̂ᵀ𝗚⁻¹∂²𝒒∂η² = 𝒒̂ᵀ𝗚⁻¹*∂²𝒒∂η²

            ∂𝒒̂ᵀ∂ξ𝗚⁻¹𝒒 =  ∂𝒒̂ᵀ∂ξ𝗚⁻¹*𝒒
            ∂𝒒̂ᵀ∂η𝗚⁻¹𝒒 =  ∂𝒒̂ᵀ∂η𝗚⁻¹*𝒒
            ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂ξ = ∂𝒒̂ᵀ∂ξ𝗚⁻¹*∂𝒒∂ξ
            ∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂ξ = ∂𝒒̂ᵀ∂η𝗚⁻¹*∂𝒒∂ξ
            ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂η = ∂𝒒̂ᵀ∂ξ𝗚⁻¹*∂𝒒∂η
            ∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂η = ∂𝒒̂ᵀ∂η𝗚⁻¹*∂𝒒∂η
            ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂²𝒒∂ξ² = ∂𝒒̂ᵀ∂ξ𝗚⁻¹*∂²𝒒∂ξ²
            ∂𝒒̂ᵀ∂η𝗚⁻¹∂²𝒒∂ξ² = ∂𝒒̂ᵀ∂η𝗚⁻¹*∂²𝒒∂ξ²
            ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂²𝒒∂ξ∂η = ∂𝒒̂ᵀ∂ξ𝗚⁻¹*∂²𝒒∂ξ∂η
            ∂𝒒̂ᵀ∂η𝗚⁻¹∂²𝒒∂ξ∂η = ∂𝒒̂ᵀ∂η𝗚⁻¹*∂²𝒒∂ξ∂η
            ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂²𝒒∂η² = ∂𝒒̂ᵀ∂ξ𝗚⁻¹*∂²𝒒∂η²
            ∂𝒒̂ᵀ∂η𝗚⁻¹∂²𝒒∂η² = ∂𝒒̂ᵀ∂η𝗚⁻¹*∂²𝒒∂η²

            ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒 = -(∂𝒒̂ᵀ∂ξ𝗚⁻¹𝒒*n₁₁+∂𝒒̂ᵀ∂η𝗚⁻¹𝒒*n₂₁)/2/𝐴
            ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒 = -(∂𝒒̂ᵀ∂ξ𝗚⁻¹𝒒*n₁₂+∂𝒒̂ᵀ∂η𝗚⁻¹𝒒*n₂₂)/2/𝐴
            ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂ξ = -(∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂ξ*n₁₁+∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂ξ*n₂₁)/2/𝐴
            ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂ξ = -(∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂ξ*n₁₂+∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂ξ*n₂₂)/2/𝐴
            ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂η = -(∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂η*n₁₁+∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂η*n₂₁)/2/𝐴
            ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂η = -(∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂η*n₁₂+∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂η*n₂₂)/2/𝐴
            ∂𝒒̂ᵀ∂x𝗚⁻¹∂²𝒒∂ξ²  = -(∂𝒒̂ᵀ∂ξ𝗚⁻¹∂²𝒒∂ξ²*n₁₁+∂𝒒̂ᵀ∂η𝗚⁻¹∂²𝒒∂ξ²*n₂₁)/2/𝐴
            ∂𝒒̂ᵀ∂y𝗚⁻¹∂²𝒒∂ξ²  = -(∂𝒒̂ᵀ∂ξ𝗚⁻¹∂²𝒒∂ξ²*n₁₂+∂𝒒̂ᵀ∂η𝗚⁻¹∂²𝒒∂ξ²*n₂₂)/2/𝐴
            ∂𝒒̂ᵀ∂x𝗚⁻¹∂²𝒒∂ξ∂η = -(∂𝒒̂ᵀ∂ξ𝗚⁻¹∂²𝒒∂ξ∂η*n₁₁+∂𝒒̂ᵀ∂η𝗚⁻¹∂²𝒒∂ξ∂η*n₂₁)/2/𝐴
            ∂𝒒̂ᵀ∂y𝗚⁻¹∂²𝒒∂ξ∂η = -(∂𝒒̂ᵀ∂ξ𝗚⁻¹∂²𝒒∂ξ∂η*n₁₂+∂𝒒̂ᵀ∂η𝗚⁻¹∂²𝒒∂ξ∂η*n₂₂)/2/𝐴
            ∂𝒒̂ᵀ∂x𝗚⁻¹∂²𝒒∂η²  = -(∂𝒒̂ᵀ∂ξ𝗚⁻¹∂²𝒒∂η²*n₁₁+∂𝒒̂ᵀ∂η𝗚⁻¹∂²𝒒∂η²*n₂₁)/2/𝐴
            ∂𝒒̂ᵀ∂y𝗚⁻¹∂²𝒒∂η²  = -(∂𝒒̂ᵀ∂ξ𝗚⁻¹∂²𝒒∂η²*n₁₂+∂𝒒̂ᵀ∂η𝗚⁻¹∂²𝒒∂η²*n₂₂)/2/𝐴

            q₁₁ = 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ²*n₁₁*n₁₁ + 2*𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ∂η*n₁₁*n₂₁ + 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂η²*n₂₁*n₂₁
            ∂q₁₁∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹∂²𝒒∂ξ²*n₁₁*n₁₁ + 2*∂𝒒̂ᵀ∂x𝗚⁻¹∂²𝒒∂ξ∂η*n₁₁*n₂₁ + ∂𝒒̂ᵀ∂x𝗚⁻¹∂²𝒒∂η²*n₂₁*n₂₁
            ∂q₁₁∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹∂²𝒒∂ξ²*n₁₁*n₁₁ + 2*∂𝒒̂ᵀ∂y𝗚⁻¹∂²𝒒∂ξ∂η*n₁₁*n₂₁ + ∂𝒒̂ᵀ∂y𝗚⁻¹∂²𝒒∂η²*n₂₁*n₂₁

            q₁₂ = 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ²*n₁₁*n₁₂ + 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ∂η*(n₁₁*n₂₂+n₁₂*n₂₁) + 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂η²*n₂₁*n₂₂
            ∂q₁₂∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹∂²𝒒∂ξ²*n₁₁*n₁₂ + ∂𝒒̂ᵀ∂x𝗚⁻¹∂²𝒒∂ξ∂η*(n₁₁*n₂₂+n₁₂*n₂₁) + ∂𝒒̂ᵀ∂x𝗚⁻¹∂²𝒒∂η²*n₂₁*n₂₂
            ∂q₁₂∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹∂²𝒒∂ξ²*n₁₁*n₁₂ + ∂𝒒̂ᵀ∂y𝗚⁻¹∂²𝒒∂ξ∂η*(n₁₁*n₂₂+n₁₂*n₂₁) + ∂𝒒̂ᵀ∂y𝗚⁻¹∂²𝒒∂η²*n₂₁*n₂₂

            q₂₂ = 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ²*n₁₂*n₁₂ + 2*𝒒̂ᵀ𝗚⁻¹∂²𝒒∂ξ∂η*n₁₂*n₂₂ + 𝒒̂ᵀ𝗚⁻¹∂²𝒒∂η²*n₂₂*n₂₂
            ∂q₂₂∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹∂²𝒒∂ξ²*n₁₂*n₁₂ + 2*∂𝒒̂ᵀ∂x𝗚⁻¹∂²𝒒∂ξ∂η*n₁₂*n₂₂ + ∂𝒒̂ᵀ∂x𝗚⁻¹∂²𝒒∂η²*n₂₂*n₂₂
            ∂q₂₂∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹∂²𝒒∂ξ²*n₁₂*n₁₂ + 2*∂𝒒̂ᵀ∂y𝗚⁻¹∂²𝒒∂ξ∂η*n₁₂*n₂₂ + ∂𝒒̂ᵀ∂y𝗚⁻¹∂²𝒒∂η²*n₂₂*n₂₂

            q₁ = 𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ*n₁₁ + 𝒒̂ᵀ𝗚⁻¹∂𝒒∂η*n₂₁
            ∂q₁∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂ξ*n₁₁ + ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂η*n₂₁
            ∂q₁∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂ξ*n₁₁ + ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂η*n₂₁
            q₂ = 𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ*n₁₂ + 𝒒̂ᵀ𝗚⁻¹∂𝒒∂η*n₂₂
            ∂q₂∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂ξ*n₁₂ + ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂η*n₂₂
            ∂q₂∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂ξ*n₁₂ + ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂η*n₂₂

            q₁n₁ = 0.0;q₂n₂ = 0.0;q₁n₂ = 0.0;q₂n₁ = 0.0
            ∂q₁∂xn₁ = 0.0;∂q₂∂xn₂ = 0.0;∂q₁∂xn₂ = 0.0;∂q₂∂xn₁ = 0.0
            ∂q₁∂yn₁ = 0.0;∂q₂∂yn₂ = 0.0;∂q₁∂yn₂ = 0.0;∂q₂∂yn₁ = 0.0
            mn₁₁n₁ = 0.0;mn₁₁n₂ = 0.0;mn₁₂n₁ = 0.0;mn₁₂n₂ = 0.0;mn₂₂n₁ = 0.0;mn₂₂n₂ = 0.0
            ∂mn₁₁∂xn₁ = 0.0;∂mn₁₁∂xn₂ = 0.0;∂mn₁₂∂xn₁ = 0.0;∂mn₁₂∂xn₂ = 0.0;∂mn₂₂∂xn₁ = 0.0;∂mn₂₂∂xn₂ = 0.0
            ∂mn₁₁∂yn₁ = 0.0;∂mn₁₁∂yn₂ = 0.0;∂mn₁₂∂yn₁ = 0.0;∂mn₁₂∂yn₂ = 0.0;∂mn₂₂∂yn₁ = 0.0;∂mn₂₂∂yn₂ = 0.0
            ms₁₁ = 0.0;ms₁₂ = 0.0;ms₂₂ = 0.0
            ∂ms₁₁∂x = 0.0;∂ms₁₂∂x = 0.0;∂ms₂₂∂x = 0.0
            ∂ms₁₁∂y = 0.0;∂ms₁₂∂y = 0.0;∂ms₂₂∂y = 0.0
            Δms₁₁ = 0.0;Δms₁₂ = 0.0;Δms₂₂ = 0.0
            Δ∂ms₁₁∂x = 0.0;Δ∂ms₁₂∂x = 0.0;Δ∂ms₂₂∂x = 0.0
            Δ∂ms₁₁∂y = 0.0;Δ∂ms₁₂∂y = 0.0;Δ∂ms₂₂∂y = 0.0
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

                ∂q₁∂xn₁ += ∂q₁∂x*n₁₁
                ∂q₁∂yn₁ += ∂q₁∂y*n₁₁
                ∂q₁∂xn₂ += ∂q₁∂x*n₁₂
                ∂q₁∂yn₂ += ∂q₁∂y*n₁₂
                ∂q₂∂xn₁ += ∂q₂∂x*n₁₁
                ∂q₂∂yn₁ += ∂q₂∂y*n₁₁
                ∂q₂∂xn₂ += ∂q₂∂x*n₁₂
                ∂q₂∂yn₂ += ∂q₂∂y*n₁₂
                ∂mn₁₁∂xn₁ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₁₁*n₁₁*n₁₁/𝐿₁²
                ∂mn₁₁∂yn₁ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₁₁*n₁₁*n₁₁/𝐿₁²
                ∂mn₁₁∂xn₂ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₁₁*n₁₁*n₁₂/𝐿₁²
                ∂mn₁₁∂yn₂ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₁₁*n₁₁*n₁₂/𝐿₁²
                ∂mn₁₂∂xn₁ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₁₁*n₁₂*n₁₁/𝐿₁²
                ∂mn₁₂∂yn₁ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₁₁*n₁₂*n₁₁/𝐿₁²
                ∂mn₁₂∂xn₂ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₁₁*n₁₂*n₁₂/𝐿₁²
                ∂mn₁₂∂yn₂ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₁₁*n₁₂*n₁₂/𝐿₁²
                ∂mn₂₂∂xn₁ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₁₂*n₁₂*n₁₁/𝐿₁²
                ∂mn₂₂∂yn₁ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₁₂*n₁₂*n₁₁/𝐿₁²
                ∂mn₂₂∂xn₂ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₁₂*n₁₂*n₁₂/𝐿₁²
                ∂mn₂₂∂yn₂ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₁₂*n₁₂*n₁₂/𝐿₁²
                ∂ms₁₁∂x += (∂q₁∂x*s₁₁+∂q₂∂x*s₁₂)*n₁₁*s₁₁/𝐿₁²
                ∂ms₁₁∂y += (∂q₁∂y*s₁₁+∂q₂∂y*s₁₂)*n₁₁*s₁₁/𝐿₁²
                ∂ms₁₂∂x += (∂q₁∂x*s₁₁+∂q₂∂x*s₁₂)*0.5*(n₁₁*s₁₂+n₁₂*s₁₁)/𝐿₁²
                ∂ms₁₂∂y += (∂q₁∂y*s₁₁+∂q₂∂y*s₁₂)*0.5*(n₁₁*s₁₂+n₁₂*s₁₁)/𝐿₁²
                ∂ms₂₂∂x += (∂q₁∂x*s₁₁+∂q₂∂x*s₁₂)*n₁₂*s₁₂/𝐿₁²
                ∂ms₂₂∂y += (∂q₁∂y*s₁₁+∂q₂∂y*s₁₂)*n₁₂*s₁₂/𝐿₁²
                if ξ.η == 0.0
                    Δms₁₁ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₁₁*s₁₁/𝐿₁²-n₂₁*s₂₁/𝐿₂²)
                    Δms₂₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₁₂*s₁₂/𝐿₁²-n₂₂*s₂₂/𝐿₂²)
                    Δms₁₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*0.5*((n₁₁*s₁₂+n₁₂*s₁₁)/𝐿₁²-(n₂₁*s₂₂+n₂₂*s₂₁)/𝐿₂²)

                    Δ∂ms₁₁∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*(n₁₁*s₁₁/𝐿₁²-n₂₁*s₂₁/𝐿₂²)
                    Δ∂ms₁₁∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*(n₁₁*s₁₁/𝐿₁²-n₂₁*s₂₁/𝐿₂²)
                    Δ∂ms₂₂∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*(n₁₂*s₁₂/𝐿₁²-n₂₂*s₂₂/𝐿₂²)
                    Δ∂ms₂₂∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*(n₁₂*s₁₂/𝐿₁²-n₂₂*s₂₂/𝐿₂²)
                    Δ∂ms₁₂∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*0.5*((n₁₁*s₁₂+n₁₂*s₁₁)/𝐿₁²-(n₂₁*s₂₂+n₂₂*s₂₁)/𝐿₂²)
                    Δ∂ms₁₂∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*0.5*((n₁₁*s₁₂+n₁₂*s₁₁)/𝐿₁²-(n₂₁*s₂₂+n₂₂*s₂₁)/𝐿₂²)
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

                ∂q₁∂xn₁ += ∂q₁∂x*n₂₁
                ∂q₁∂yn₁ += ∂q₁∂y*n₂₁
                ∂q₁∂xn₂ += ∂q₁∂x*n₂₂
                ∂q₁∂yn₂ += ∂q₁∂y*n₂₂
                ∂q₂∂xn₁ += ∂q₂∂x*n₂₁
                ∂q₂∂yn₁ += ∂q₂∂y*n₂₁
                ∂q₂∂xn₂ += ∂q₂∂x*n₂₂
                ∂q₂∂yn₂ += ∂q₂∂y*n₂₂
                ∂mn₁₁∂xn₁ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₂₁*n₂₁*n₂₁/𝐿₂²
                ∂mn₁₁∂yn₁ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₂₁*n₂₁*n₂₁/𝐿₂²
                ∂mn₁₁∂xn₂ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₂₁*n₂₁*n₂₂/𝐿₂²
                ∂mn₁₁∂yn₂ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₂₁*n₂₁*n₂₂/𝐿₂²
                ∂mn₁₂∂xn₁ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₂₁*n₂₂*n₂₁/𝐿₂²
                ∂mn₁₂∂yn₁ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₂₁*n₂₂*n₂₁/𝐿₂²
                ∂mn₁₂∂xn₂ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₂₁*n₂₂*n₂₂/𝐿₂²
                ∂mn₁₂∂yn₂ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₂₁*n₂₂*n₂₂/𝐿₂²
                ∂mn₂₂∂xn₁ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₂₂*n₂₂*n₂₁/𝐿₂²
                ∂mn₂₂∂yn₁ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₂₂*n₂₂*n₂₁/𝐿₂²
                ∂mn₂₂∂xn₂ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₂₂*n₂₂*n₂₂/𝐿₂²
                ∂mn₂₂∂yn₂ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₂₂*n₂₂*n₂₂/𝐿₂²
                ∂ms₁₁∂x += (∂q₁∂x*s₂₁+∂q₂∂x*s₂₂)*n₂₁*s₂₁/𝐿₂²
                ∂ms₁₁∂y += (∂q₁∂y*s₂₁+∂q₂∂y*s₂₂)*n₂₁*s₂₁/𝐿₂²
                ∂ms₁₂∂x += (∂q₁∂x*s₂₁+∂q₂∂x*s₂₂)*0.5*(n₂₁*s₂₂+n₂₂*s₂₁)/𝐿₂²
                ∂ms₁₂∂y += (∂q₁∂y*s₂₁+∂q₂∂y*s₂₂)*0.5*(n₂₁*s₂₂+n₂₂*s₂₁)/𝐿₂²
                ∂ms₂₂∂x += (∂q₁∂x*s₂₁+∂q₂∂x*s₂₂)*n₂₂*s₂₂/𝐿₂²
                ∂ms₂₂∂y += (∂q₁∂y*s₂₁+∂q₂∂y*s₂₂)*n₂₂*s₂₂/𝐿₂²

                if ξ.ξ+ξ.η ≈ 1.0
                    Δms₁₁ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₂₁*s₂₁/𝐿₂²-n₃₁*s₃₁/𝐿₃²)
                    Δms₂₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₂₂*s₂₂/𝐿₂²-n₃₂*s₃₂/𝐿₃²)
                    Δms₁₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*0.5*((n₂₁*s₂₂+n₂₂*s₂₁)/𝐿₂²-(n₃₁*s₃₂+n₃₂*s₃₁)/𝐿₃²)

                    Δ∂ms₁₁∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*(n₂₁*s₂₁/𝐿₂²-n₃₁*s₃₁/𝐿₃²)
                    Δ∂ms₁₁∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*(n₂₁*s₂₁/𝐿₂²-n₃₁*s₃₁/𝐿₃²)
                    Δ∂ms₂₂∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*(n₂₂*s₂₂/𝐿₂²-n₃₂*s₃₂/𝐿₃²)
                    Δ∂ms₂₂∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*(n₂₂*s₂₂/𝐿₂²-n₃₂*s₃₂/𝐿₃²)
                    Δ∂ms₁₂∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*0.5*((n₂₁*s₂₂+n₂₂*s₂₁)/𝐿₂²-(n₃₁*s₃₂+n₃₂*s₃₁)/𝐿₃²)
                    Δ∂ms₁₂∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*0.5*((n₂₁*s₂₂+n₂₂*s₂₁)/𝐿₂²-(n₃₁*s₃₂+n₃₂*s₃₁)/𝐿₃²)
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

                ∂q₁∂xn₁ += ∂q₁∂x*n₃₁
                ∂q₁∂yn₁ += ∂q₁∂y*n₃₁
                ∂q₁∂xn₂ += ∂q₁∂x*n₃₂
                ∂q₁∂yn₂ += ∂q₁∂y*n₃₂
                ∂q₂∂xn₁ += ∂q₂∂x*n₃₁
                ∂q₂∂yn₁ += ∂q₂∂y*n₃₁
                ∂q₂∂xn₂ += ∂q₂∂x*n₃₂
                ∂q₂∂yn₂ += ∂q₂∂y*n₃₂
                ∂mn₁₁∂xn₁ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₃₁*n₃₁*n₃₁/𝐿₃²
                ∂mn₁₁∂yn₁ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₃₁*n₃₁*n₃₁/𝐿₃²
                ∂mn₁₁∂xn₂ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₃₁*n₃₁*n₃₂/𝐿₃²
                ∂mn₁₁∂yn₂ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₃₁*n₃₁*n₃₂/𝐿₃²
                ∂mn₁₂∂xn₁ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₃₁*n₃₂*n₃₁/𝐿₃²
                ∂mn₁₂∂yn₁ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₃₁*n₃₂*n₃₁/𝐿₃²
                ∂mn₁₂∂xn₂ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₃₁*n₃₂*n₃₂/𝐿₃²
                ∂mn₁₂∂yn₂ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₃₁*n₃₂*n₃₂/𝐿₃²
                ∂mn₂₂∂xn₁ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₃₂*n₃₂*n₃₁/𝐿₃²
                ∂mn₂₂∂yn₁ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₃₂*n₃₂*n₃₁/𝐿₃²
                ∂mn₂₂∂xn₂ += ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₃₂*n₃₂*n₃₂/𝐿₃²
                ∂mn₂₂∂yn₂ += ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₃₂*n₃₂*n₃₂/𝐿₃²
                ∂ms₁₁∂x += (∂q₁∂x*s₃₁+∂q₂∂x*s₃₂)*n₃₁*s₃₁/𝐿₃²
                ∂ms₁₁∂y += (∂q₁∂y*s₃₁+∂q₂∂y*s₃₂)*n₃₁*s₃₁/𝐿₃²
                ∂ms₁₂∂x += (∂q₁∂x*s₃₁+∂q₂∂x*s₃₂)*0.5*(n₃₁*s₃₂+n₃₂*s₃₁)/𝐿₃²
                ∂ms₁₂∂y += (∂q₁∂y*s₃₁+∂q₂∂y*s₃₂)*0.5*(n₃₁*s₃₂+n₃₂*s₃₁)/𝐿₃²
                ∂ms₂₂∂x += (∂q₁∂x*s₃₁+∂q₂∂x*s₃₂)*n₃₂*s₃₂/𝐿₃²
                ∂ms₂₂∂y += (∂q₁∂y*s₃₁+∂q₂∂y*s₃₂)*n₃₂*s₃₂/𝐿₃²
                if ξ.ξ == 0.0
                    Δms₁₁ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₃₁*s₃₁/𝐿₃²-n₁₁*s₁₁/𝐿₁²)
                    Δms₂₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*(n₃₂*s₃₂/𝐿₃²-n₁₂*s₁₂/𝐿₁²)
                    Δms₁₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*0.5*((n₃₁*s₃₂+n₃₂*s₃₁)/𝐿₃²-(n₁₁*s₁₂+n₁₂*s₁₁)/𝐿₁²)

                    Δ∂ms₁₁∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*(n₃₁*s₃₁/𝐿₃²-n₁₁*s₁₁/𝐿₁²)
                    Δ∂ms₁₁∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*(n₃₁*s₃₁/𝐿₃²-n₁₁*s₁₁/𝐿₁²)
                    Δ∂ms₂₂∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*(n₃₂*s₃₂/𝐿₃²-n₁₂*s₁₂/𝐿₁²)
                    Δ∂ms₂₂∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*(n₃₂*s₃₂/𝐿₃²-n₁₂*s₁₂/𝐿₁²)
                    Δ∂ms₁₂∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*0.5*((n₃₁*s₃₂+n₃₂*s₃₁)/𝐿₃²-(n₁₁*s₁₂+n₁₂*s₁₁)/𝐿₁²)
                    Δ∂ms₁₂∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*0.5*((n₃₁*s₃₂+n₃₂*s₃₁)/𝐿₃²-(n₁₁*s₁₂+n₁₂*s₁₁)/𝐿₁²)
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

            ∂W₁₁₁∂x = ∂mn₁₁∂xn₁*wᵇ
            ∂W₁₁₁∂y = ∂mn₁₁∂yn₁*wᵇ
            ∂W₁₁₂∂x = ∂mn₁₁∂xn₂*wᵇ
            ∂W₁₁₂∂y = ∂mn₁₁∂yn₂*wᵇ
            ∂W₁₂₁∂x = ∂mn₁₂∂xn₁*wᵇ
            ∂W₁₂₁∂y = ∂mn₁₂∂yn₁*wᵇ
            ∂W₁₂₂∂x = ∂mn₁₂∂xn₂*wᵇ
            ∂W₁₂₂∂y = ∂mn₁₂∂yn₂*wᵇ
            ∂W₂₂₁∂x = ∂mn₂₂∂xn₁*wᵇ
            ∂W₂₂₁∂y = ∂mn₂₂∂yn₁*wᵇ
            ∂W₂₂₂∂x = ∂mn₂₂∂xn₂*wᵇ
            ∂W₂₂₂∂y = ∂mn₂₂∂yn₂*wᵇ
            ∂W₁₁∂x = (∂q₁₁∂x*w + 2*(∂q₁∂xn₁+∂ms₁₁∂x)*wᵇ)/4/𝐴 + Δ∂ms₁₁∂x
            ∂W₁₁∂y = (∂q₁₁∂y*w + 2*(∂q₁∂yn₁+∂ms₁₁∂y)*wᵇ)/4/𝐴 + Δ∂ms₁₁∂y
            ∂W₁₂∂x = (∂q₁₂∂x*w + (∂q₁∂xn₂+∂q₂∂xn₁+2*∂ms₁₂∂x)*wᵇ)/4/𝐴 + Δ∂ms₁₂∂x
            ∂W₁₂∂y = (∂q₁₂∂y*w + (∂q₁∂yn₂+∂q₂∂yn₁+2*∂ms₁₂∂y)*wᵇ)/4/𝐴 + Δ∂ms₁₂∂y
            ∂W₂₂∂x = (∂q₂₂∂x*w + 2*(∂q₂∂xn₂+∂ms₂₂∂x)*wᵇ)/4/𝐴 + Δ∂ms₂₂∂x
            ∂W₂₂∂y = (∂q₂₂∂y*w + 2*(∂q₂∂yn₂+∂ms₂₂∂y)*wᵇ)/4/𝐴 + Δ∂ms₂₂∂y
            for i in 1:length(𝓒)
                ∂²𝝭∂x²[i] += 𝝭[i]*W₁₁ + ∂𝝭∂x[i]*W₁₁₁ + ∂𝝭∂y[i]*W₁₁₂
                ∂²𝝭∂x∂y[i] += 𝝭[i]*W₁₂ + ∂𝝭∂x[i]*W₁₂₁ + ∂𝝭∂y[i]*W₁₂₂
                ∂²𝝭∂y²[i] += 𝝭[i]*W₂₂ + ∂𝝭∂x[i]*W₂₂₁ + ∂𝝭∂y[i]*W₂₂₂
                ∂∂²𝝭∂x²∂x[i] += 𝝭[i]*∂W₁₁∂x + ∂𝝭∂x[i]*∂W₁₁₁∂x + ∂𝝭∂y[i]*∂W₁₁₂∂x
                ∂∂²𝝭∂x²∂y[i] += 𝝭[i]*∂W₁₁∂y + ∂𝝭∂x[i]*∂W₁₁₁∂y + ∂𝝭∂y[i]*∂W₁₁₂∂y
                ∂∂²𝝭∂x∂y∂x[i] += 𝝭[i]*∂W₁₂∂x + ∂𝝭∂x[i]*∂W₁₂₁∂x + ∂𝝭∂y[i]*∂W₁₂₂∂x
                ∂∂²𝝭∂x∂y∂y[i] += 𝝭[i]*∂W₁₂∂y + ∂𝝭∂x[i]*∂W₁₂₁∂y + ∂𝝭∂y[i]*∂W₁₂₂∂y
                ∂∂²𝝭∂y²∂x[i] += 𝝭[i]*∂W₂₂∂x + ∂𝝭∂x[i]*∂W₂₂₁∂x + ∂𝝭∂y[i]*∂W₂₂₂∂x
                ∂∂²𝝭∂y²∂y[i] += 𝝭[i]*∂W₂₂∂y + ∂𝝭∂x[i]*∂W₂₂₁∂y + ∂𝝭∂y[i]*∂W₂₂₂∂y
            end
        end
    end
end

function set∇̄𝝭!(ap::ReproducingKernel{𝒑,𝑠,𝜙,:Seg2}) where {𝒑,𝑠,𝜙}
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ̂ in 𝓖
        𝒒̂ = get𝒑₁(ap,ξ̂)
        𝗚⁻¹ = cal𝗚!(ap)
        𝒒̂ᵀ𝗚⁻¹ = 𝒒̂*𝗚⁻¹
        ∂𝝭∂x = ξ̂[:∂𝝭∂x_]
        fill!(∂𝝭∂x,0.0)
        for ξ in ap.𝓖
            w = ξ.w
            n = ξ.n₁
            𝝭 = ξ[:𝝭]
            𝒒 = get𝒑₁(ap,ξ)
            W₁ = 𝒒̂ᵀ𝗚⁻¹*𝒒*n*w
            for i in 1:length(𝓒)
                ∂𝝭∂x[i] += 𝝭[i]*W₁
            end
        end
    end
end

function set∇̄𝝭!(ap::ReproducingKernel{𝒑,𝑠,𝜙,:Tri3}) where {𝒑,𝑠,𝜙}
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ̂ in 𝓖
        𝒒̂ = get𝒑₁(ap,ξ̂)
        𝗚⁻¹ = cal𝗚!(ap)
        𝒒̂ᵀ𝗚⁻¹ = 𝒒̂*𝗚⁻¹
        ∂𝝭∂x = ξ̂[:∂𝝭∂x_]
        ∂𝝭∂y = ξ̂[:∂𝝭∂y_]
        for i in 1:length(𝓒)
            ∂𝝭∂x[i] = 0.0
            ∂𝝭∂y[i] = 0.0
        end
        for ξ in ap.𝓖
            w = ξ.w
            n₁ = ξ.n₁
            n₂ = ξ.n₂
            𝝭 = ξ[:𝝭]
            𝒒 = get𝒑₁(ap,ξ)
            𝒒̂ᵀ𝗚⁻¹𝒒 = 𝒒̂ᵀ𝗚⁻¹*𝒒
            W₁ = 𝒒̂ᵀ𝗚⁻¹𝒒*n₁*w
            W₂ = 𝒒̂ᵀ𝗚⁻¹𝒒*n₂*w
            for i in 1:length(𝓒)
                ∂𝝭∂x[i] += 𝝭[i]*W₁
                ∂𝝭∂y[i] += 𝝭[i]*W₂
            end
        end
    end
end

function set∇∇̄²𝝭!(ap::ReproducingKernel{𝒑,𝑠,𝜙,:Tri3};Γᵍ::Union{ReproducingKernel{𝒑,𝑠,𝜙,:Tri3},Nothing}=nothing,Γᶿ::Union{ReproducingKernel{𝒑,𝑠,𝜙,:Tri3},Nothing}=nothing,Γᴾ::Vector{ReproducingKernel{𝒑,𝑠,𝜙,:Tri3}}=[]) where {𝒑,𝑠,𝜙}
    x₁ = ap.𝓒[1].x;y₁ = ap.𝓒[1].y
    x₂ = ap.𝓒[2].x;y₂ = ap.𝓒[2].y
    x₃ = ap.𝓒[3].x;y₃ = ap.𝓒[3].y
    𝐴 = get𝐴(ap)
    n₁₁ = y₃-y₂;n₂₁ = y₁-y₃;n₃₁ = y₂-y₁
    n₁₂ = x₂-x₃;n₂₂ = x₃-x₁;n₃₂ = x₁-x₂
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ̂ in 𝓖
        𝒒̂,∂𝒒̂∂ξ,∂𝒒̂∂η = get∇𝒑₂(ap,ξ̂)
        𝗚⁻¹,𝗚⁻¹∂ξ,𝗚⁻¹∂η = cal∇𝗚₂!(ap)
        𝒒̂ᵀ𝗚⁻¹ = 𝒒̂*𝗚⁻¹
        ∂𝒒̂ᵀ∂ξ𝗚⁻¹ = ∂𝒒̂∂ξ*𝗚⁻¹∂ξ
        ∂𝒒̂ᵀ∂η𝗚⁻¹ = ∂𝒒̂∂η*𝗚⁻¹∂η

        ∂²𝝭∂x² = ξ̂[:∂²𝝭∂x²_]
        ∂²𝝭∂x∂y = ξ̂[:∂²𝝭∂x∂y_]
        ∂²𝝭∂y² = ξ̂[:∂²𝝭∂y²_]
        ∂∂²𝝭∂x²∂x = ξ̂[:∂∂²𝝭∂x²∂x_]
        ∂∂²𝝭∂x²∂y = ξ̂[:∂∂²𝝭∂x²∂y_]
        ∂∂²𝝭∂x∂y∂x = ξ̂[:∂∂²𝝭∂x∂y∂x_]
        ∂∂²𝝭∂x∂y∂y = ξ̂[:∂∂²𝝭∂x∂y∂y_]
        ∂∂²𝝭∂y²∂x = ξ̂[:∂∂²𝝭∂y²∂x_]
        ∂∂²𝝭∂y²∂y = ξ̂[:∂∂²𝝭∂y²∂y_]
        for i in 1:length(𝓒)
            ∂²𝝭∂x²[i] = 0.0
            ∂²𝝭∂x∂y[i] = 0.0
            ∂²𝝭∂y²[i] = 0.0
            ∂∂²𝝭∂x²∂x[i] = 0.0
            ∂∂²𝝭∂x²∂y[i] = 0.0
            ∂∂²𝝭∂x∂y∂x[i] = 0.0
            ∂∂²𝝭∂x∂y∂y[i] = 0.0
            ∂∂²𝝭∂y²∂x[i] = 0.0
            ∂∂²𝝭∂y²∂y[i] = 0.0
        end
        if Γᵍ ≠ nothing
            for ξ in Γᵍ.𝓖
                𝑤 = ξ.𝑤
                n₁ = ξ.n₁
                n₂ = ξ.n₂
                s₁ = -n₂
                s₂ = n₁
                𝝭 = ξ[:𝝭]
                𝒒, ∂𝒒∂ξ, ∂𝒒∂η = get∇²𝒑₂(Γᵍ,ξ)

                𝒒̂ᵀ𝗚⁻¹𝒒 =  𝒒̂ᵀ𝗚⁻¹*𝒒
                𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ = 𝒒̂ᵀ𝗚⁻¹*∂𝒒∂ξ
                𝒒̂ᵀ𝗚⁻¹∂𝒒∂η = 𝒒̂ᵀ𝗚⁻¹*∂𝒒∂η

                ∂𝒒̂ᵀ∂ξ𝗚⁻¹𝒒 =  ∂𝒒̂ᵀ∂ξ𝗚⁻¹*𝒒
                ∂𝒒̂ᵀ∂η𝗚⁻¹𝒒 =  ∂𝒒̂ᵀ∂η𝗚⁻¹*𝒒
                ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂ξ = ∂𝒒̂ᵀ∂ξ𝗚⁻¹*∂𝒒∂ξ
                ∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂ξ = ∂𝒒̂ᵀ∂η𝗚⁻¹*∂𝒒∂ξ
                ∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂η = ∂𝒒̂ᵀ∂ξ𝗚⁻¹*∂𝒒∂η
                ∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂η = ∂𝒒̂ᵀ∂η𝗚⁻¹*∂𝒒∂η

                ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂ξ = -(∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂ξ*n₁₁+∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂ξ*n₂₁)/2/𝐴
                ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂ξ = -(∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂ξ*n₁₂+∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂ξ*n₂₂)/2/𝐴
                ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂η = -(∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂η*n₁₁+∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂η*n₂₁)/2/𝐴
                ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂η = -(∂𝒒̂ᵀ∂ξ𝗚⁻¹∂𝒒∂η*n₁₂+∂𝒒̂ᵀ∂η𝗚⁻¹∂𝒒∂η*n₂₂)/2/𝐴

                q₁ = 𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ*n₁₁ + 𝒒̂ᵀ𝗚⁻¹∂𝒒∂η*n₂₁
                ∂q₁∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂ξ*n₁₁ + ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂η*n₂₁
                ∂q₁∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂ξ*n₁₁ + ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂η*n₂₁
                q₂ = 𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ*n₁₂ + 𝒒̂ᵀ𝗚⁻¹∂𝒒∂η*n₂₂
                ∂q₂∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂ξ*n₁₂ + ∂𝒒̂ᵀ∂x𝗚⁻¹∂𝒒∂η*n₂₂
                ∂q₂∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂ξ*n₁₂ + ∂𝒒̂ᵀ∂y𝗚⁻¹∂𝒒∂η*n₂₂

                q₁n₁ = q₁*n₁
                q₁n₂ = q₁*n₂
                q₂n₁ = q₂*n₁
                q₂n₂ = q₂*n₂
                ms₁₁ = (q₁*s₁+q₂*s₂)*n₁*s₁
                ms₁₂ = (q₁*s₁+q₂*s₂)*0.5*(n₁*s₂+n₂*s₁)
                ms₂₂ = (q₁*s₁+q₂*s₂)*n₂*s₂

                ∂q₁∂xn₁ = ∂q₁∂x*n₁
                ∂q₁∂yn₁ = ∂q₁∂y*n₁
                ∂q₁∂xn₂ = ∂q₁∂x*n₂
                ∂q₁∂yn₂ = ∂q₁∂y*n₂
                ∂q₂∂xn₁ = ∂q₂∂x*n₁
                ∂q₂∂yn₁ = ∂q₂∂y*n₁
                ∂q₂∂xn₂ = ∂q₂∂x*n₂
                ∂q₂∂yn₂ = ∂q₂∂y*n₂
                ∂ms₁₁∂x = (∂q₁∂x*s₁+∂q₂∂x*s₂)*n₁*s₁
                ∂ms₁₁∂y = (∂q₁∂y*s₁+∂q₂∂y*s₂)*n₁*s₁
                ∂ms₁₂∂x = (∂q₁∂x*s₁+∂q₂∂x*s₂)*0.5*(n₁*s₂+n₂*s₁)
                ∂ms₁₂∂y = (∂q₁∂y*s₁+∂q₂∂y*s₂)*0.5*(n₁*s₂+n₂*s₁)
                ∂ms₂₂∂x = (∂q₁∂x*s₁+∂q₂∂x*s₂)*n₂*s₂
                ∂ms₂₂∂y = (∂q₁∂y*s₁+∂q₂∂y*s₂)*n₂*s₂

                W₁₁ = (q₁n₁+ms₁₁)*𝑤/2/𝐴
                W₁₂ = (q₁n₂+q₂n₁+2*ms₁₂)*𝑤/4/𝐴
                W₂₂ = (q₂n₂+ms₂₂)*𝑤/2/𝐴
                ∂W₁₁∂x = (∂q₁∂xn₁+∂ms₁₁∂x)*𝑤/2/𝐴
                ∂W₁₁∂y = (∂q₁∂yn₁+∂ms₁₁∂y)*𝑤/2/𝐴
                ∂W₁₂∂x = (∂q₁∂xn₂+∂q₂∂xn₁+2*∂ms₁₂∂x)*𝑤/4/𝐴
                ∂W₁₂∂y = (∂q₁∂yn₂+∂q₂∂yn₁+2*∂ms₁₂∂y)*𝑤/4/𝐴
                ∂W₂₂∂x = (∂q₂∂xn₂+∂ms₂₂∂x)*𝑤/2/𝐴
                ∂W₂₂∂y = (∂q₂∂yn₂+∂ms₂₂∂y)*𝑤/2/𝐴
                for i in 1:length(𝓒)
                    ∂²𝝭∂x²[i] += 𝝭[i]*W₁₁
                    ∂²𝝭∂x∂y[i] += 𝝭[i]*W₁₂
                    ∂²𝝭∂y²[i] += 𝝭[i]*W₂₂
                    ∂∂²𝝭∂x²∂x[i] += 𝝭[i]*∂W₁₁∂x
                    ∂∂²𝝭∂x²∂y[i] += 𝝭[i]*∂W₁₁∂y
                    ∂∂²𝝭∂x∂y∂x[i] += 𝝭[i]*∂W₁₂∂x
                    ∂∂²𝝭∂x∂y∂y[i] += 𝝭[i]*∂W₁₂∂y
                    ∂∂²𝝭∂y²∂x[i] += 𝝭[i]*∂W₂₂∂x
                    ∂∂²𝝭∂y²∂y[i] += 𝝭[i]*∂W₂₂∂y
                end
            end
        end

        if Γᶿ ≠ nothing
            for ξ in Γᶿ.𝓖
                𝑤 = ξ.𝑤
                n₁ = ξ.n₁
                n₂ = ξ.n₂
                s₁ = ξ.s₁
                s₂ = ξ.s₂
                ∂𝝭∂x = ξ[:∂𝝭∂x]
                ∂𝝭∂y = ξ[:∂𝝭∂y]
                𝒒 = get𝒑₂(Γᶿ,ξ)

                𝒒̂ᵀ𝗚⁻¹𝒒 =  𝒒̂ᵀ𝗚⁻¹*𝒒
                ∂𝒒̂ᵀ∂ξ𝗚⁻¹𝒒 =  ∂𝒒̂ᵀ∂ξ𝗚⁻¹*𝒒
                ∂𝒒̂ᵀ∂η𝗚⁻¹𝒒 =  ∂𝒒̂ᵀ∂η𝗚⁻¹*𝒒

                ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒 = -(∂𝒒̂ᵀ∂ξ𝗚⁻¹𝒒*n₁₁+∂𝒒̂ᵀ∂η𝗚⁻¹𝒒*n₂₁)/2/𝐴
                ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒 = -(∂𝒒̂ᵀ∂ξ𝗚⁻¹𝒒*n₁₂+∂𝒒̂ᵀ∂η𝗚⁻¹𝒒*n₂₂)/2/𝐴

                mn₁₁n₁ = 𝒒̂ᵀ𝗚⁻¹𝒒*n₁*n₁*n₁
                mn₁₁n₂ = 𝒒̂ᵀ𝗚⁻¹𝒒*n₁*n₁*n₂
                mn₁₂n₁ = 𝒒̂ᵀ𝗚⁻¹𝒒*n₁*n₂*n₁
                mn₁₂n₂ = 𝒒̂ᵀ𝗚⁻¹𝒒*n₁*n₂*n₂
                mn₂₂n₁ = 𝒒̂ᵀ𝗚⁻¹𝒒*n₂*n₂*n₁
                mn₂₂n₂ = 𝒒̂ᵀ𝗚⁻¹𝒒*n₂*n₂*n₂

                ∂mn₁₁∂xn₁ = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₁*n₁*n₁
                ∂mn₁₁∂yn₁ = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₁*n₁*n₁
                ∂mn₁₁∂xn₂ = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₁*n₁*n₂
                ∂mn₁₁∂yn₂ = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₁*n₁*n₂
                ∂mn₁₂∂xn₁ = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₁*n₂*n₁
                ∂mn₁₂∂yn₁ = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₁*n₂*n₁
                ∂mn₁₂∂xn₂ = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₁*n₂*n₂
                ∂mn₁₂∂yn₂ = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₁*n₂*n₂
                ∂mn₂₂∂xn₁ = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₂*n₂*n₁
                ∂mn₂₂∂yn₁ = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₂*n₂*n₁
                ∂mn₂₂∂xn₂ = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*n₂*n₂*n₂
                ∂mn₂₂∂yn₂ = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*n₂*n₂*n₂

                W₁₁₁ = mn₁₁n₁*𝑤
                W₁₁₂ = mn₁₁n₂*𝑤
                W₁₂₁ = mn₁₂n₁*𝑤
                W₁₂₂ = mn₁₂n₂*𝑤
                W₂₂₁ = mn₂₂n₁*𝑤
                W₂₂₂ = mn₂₂n₂*𝑤

                ∂W₁₁₁∂x = ∂mn₁₁∂xn₁*𝑤
                ∂W₁₁₁∂y = ∂mn₁₁∂yn₁*𝑤
                ∂W₁₁₂∂x = ∂mn₁₁∂xn₂*𝑤
                ∂W₁₁₂∂y = ∂mn₁₁∂yn₂*𝑤
                ∂W₁₂₁∂x = ∂mn₁₂∂xn₁*𝑤
                ∂W₁₂₁∂y = ∂mn₁₂∂yn₁*𝑤
                ∂W₁₂₂∂x = ∂mn₁₂∂xn₂*𝑤
                ∂W₁₂₂∂y = ∂mn₁₂∂yn₂*𝑤
                ∂W₂₂₁∂x = ∂mn₂₂∂xn₁*𝑤
                ∂W₂₂₁∂y = ∂mn₂₂∂yn₁*𝑤
                ∂W₂₂₂∂x = ∂mn₂₂∂xn₂*𝑤
                ∂W₂₂₂∂y = ∂mn₂₂∂yn₂*𝑤
                for i in 1:length(𝓒)
                    ∂²𝝭∂x²[i] += ∂𝝭∂x[i]*W₁₁₁ + ∂𝝭∂y[i]*W₁₁₂
                    ∂²𝝭∂x∂y[i] += ∂𝝭∂x[i]*W₁₂₁ + ∂𝝭∂y[i]*W₁₂₂
                    ∂²𝝭∂y²[i] += ∂𝝭∂x[i]*W₂₂₁ + ∂𝝭∂y[i]*W₂₂₂
                    ∂∂²𝝭∂x²∂x[i] += ∂𝝭∂x[i]*∂W₁₁₁∂x + ∂𝝭∂y[i]*∂W₁₁₂∂x
                    ∂∂²𝝭∂x²∂y[i] += ∂𝝭∂x[i]*∂W₁₁₁∂y + ∂𝝭∂y[i]*∂W₁₁₂∂y
                    ∂∂²𝝭∂x∂y∂x[i] += ∂𝝭∂x[i]*∂W₁₂₁∂x + ∂𝝭∂y[i]*∂W₁₂₂∂x
                    ∂∂²𝝭∂x∂y∂y[i] += ∂𝝭∂x[i]*∂W₁₂₁∂y + ∂𝝭∂y[i]*∂W₁₂₂∂y
                    ∂∂²𝝭∂y²∂x[i] += ∂𝝭∂x[i]*∂W₂₂₁∂x + ∂𝝭∂y[i]*∂W₂₂₂∂x
                    ∂∂²𝝭∂y²∂y[i] += ∂𝝭∂x[i]*∂W₂₂₁∂y + ∂𝝭∂y[i]*∂W₂₂₂∂y
                end
            end
        end

        for a in Γᴾ
            if ap∩a ≠ nothing
                ξ = a.𝓖[1]
                Δn₁s₁ = ξ.Δn₁s₁
                Δn₁s₂n₂s₁ = ξ.Δn₁s₂n₂s₁
                Δn₂s₂ = ξ.Δn₂s₂

                𝝭 = ξ[:𝝭]
                𝒒 = get𝒑₂(a,ξ)

                𝒒̂ᵀ𝗚⁻¹𝒒 =  𝒒̂ᵀ𝗚⁻¹*𝒒
                ∂𝒒̂ᵀ∂ξ𝗚⁻¹𝒒 =  ∂𝒒̂ᵀ∂ξ𝗚⁻¹*𝒒
                ∂𝒒̂ᵀ∂η𝗚⁻¹𝒒 =  ∂𝒒̂ᵀ∂η𝗚⁻¹*𝒒

                ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒 = -(∂𝒒̂ᵀ∂ξ𝗚⁻¹𝒒*n₁₁+∂𝒒̂ᵀ∂η𝗚⁻¹𝒒*n₂₁)/2/𝐴
                ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒 = -(∂𝒒̂ᵀ∂ξ𝗚⁻¹𝒒*n₁₂+∂𝒒̂ᵀ∂η𝗚⁻¹𝒒*n₂₂)/2/𝐴

                Δms₁₁ =  𝒒̂ᵀ𝗚⁻¹𝒒*Δn₁s₁
                Δms₂₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*Δn₂s₂
                Δms₁₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*Δn₁s₂n₂s₁/2

                Δ∂ms₁₁∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*Δn₁s₁
                Δ∂ms₁₁∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*Δn₁s₁
                Δ∂ms₂₂∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*Δn₂s₂
                Δ∂ms₂₂∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*Δn₂s₂
                Δ∂ms₁₂∂x = ∂𝒒̂ᵀ∂x𝗚⁻¹𝒒*Δn₁s₂n₂s₁/2
                Δ∂ms₁₂∂y = ∂𝒒̂ᵀ∂y𝗚⁻¹𝒒*Δn₁s₂n₂s₁/2

                for i in 1:length(𝓒)
                    ∂²𝝭∂x²[i] += 𝝭[i]*Δms₁₁
                    ∂²𝝭∂x∂y[i] += 𝝭[i]*Δms₁₂
                    ∂²𝝭∂y²[i] += 𝝭[i]*Δms₂₂
                    ∂∂²𝝭∂x²∂x[i] += 𝝭[i]*Δ∂ms₁₁∂x
                    ∂∂²𝝭∂x²∂y[i] += 𝝭[i]*Δ∂ms₁₁∂y
                    ∂∂²𝝭∂x∂y∂x[i] += 𝝭[i]*Δ∂ms₁₂∂x
                    ∂∂²𝝭∂x∂y∂y[i] += 𝝭[i]*Δ∂ms₁₂∂y
                    ∂∂²𝝭∂y²∂x[i] += 𝝭[i]*Δ∂ms₂₂∂x
                    ∂∂²𝝭∂y²∂y[i] += 𝝭[i]*Δ∂ms₂₂∂y
                end
            end
        end
    end
end

function set∇̄²𝝭!(ap::ReproducingKernel{𝒑,𝑠,𝜙,:Tri3};Γᵍ::Vector{ReproducingKernel{𝒑,𝑠,𝜙,:Tri3}}=[],Γᶿ::Vector{ReproducingKernel{𝒑,𝑠,𝜙,:Tri3}}=[],Γᴾ::Vector{ReproducingKernel{𝒑,𝑠,𝜙,:Tri3}}=[]) where {𝒑,𝑠,𝜙}
    x₁ = ap.𝓒[1].x;y₁ = ap.𝓒[1].y
    x₂ = ap.𝓒[2].x;y₂ = ap.𝓒[2].y
    x₃ = ap.𝓒[3].x;y₃ = ap.𝓒[3].y
    𝐴 = get𝐴(ap)
    n₁₁ = y₃-y₂;n₂₁ = y₁-y₃;n₃₁ = y₂-y₁
    n₁₂ = x₂-x₃;n₂₂ = x₃-x₁;n₃₂ = x₁-x₂
    s₁₁ = -n₁₂;s₂₁ = -n₂₂;s₃₁ = -n₃₂
    s₁₂ =  n₁₁;s₂₂ =  n₂₁;s₃₂ =  n₃₁
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ̂ in 𝓖
        𝒒̂ = get𝒑₂(ap,ξ̂)
        𝗚⁻¹ = cal𝗚₂!(ap)
        𝒒̂ᵀ𝗚⁻¹ = 𝒒̂*𝗚⁻¹

        ∂²𝝭∂x² = ξ̂[:∂²𝝭∂x²_]
        ∂²𝝭∂x∂y = ξ̂[:∂²𝝭∂x∂y_]
        ∂²𝝭∂y² = ξ̂[:∂²𝝭∂y²_]
        for i in 1:length(𝓒)
            ∂²𝝭∂x²[i] = 0.0
            ∂²𝝭∂x∂y[i] = 0.0
            ∂²𝝭∂y²[i] = 0.0
        end
        for a in Γᵍ
            if ap∩a ≠ nothing
                for ξ in a.𝓖
                    𝑤 = ξ.𝑤
                    n₁ = ξ.n₁
                    n₂ = ξ.n₂
                    s₁ = -n₂
                    s₂ = n₁
                    𝝭 = ξ[:𝝭]
                    𝒒, ∂𝒒∂ξ, ∂𝒒∂η = get∇²𝒑₂(a,ξ)

                    𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ = 𝒒̂ᵀ𝗚⁻¹*∂𝒒∂ξ
                    𝒒̂ᵀ𝗚⁻¹∂𝒒∂η = 𝒒̂ᵀ𝗚⁻¹*∂𝒒∂η

                    q₁ = 𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ*n₁₁ + 𝒒̂ᵀ𝗚⁻¹∂𝒒∂η*n₂₁
                    q₂ = 𝒒̂ᵀ𝗚⁻¹∂𝒒∂ξ*n₁₂ + 𝒒̂ᵀ𝗚⁻¹∂𝒒∂η*n₂₂

                    q₁n₁ = q₁*n₁
                    q₁n₂ = q₁*n₂
                    q₂n₁ = q₂*n₁
                    q₂n₂ = q₂*n₂
                    ms₁₁ = (q₁*s₁+q₂*s₂)*n₁*s₁
                    ms₁₂ = (q₁*s₁+q₂*s₂)*0.5*(n₁*s₂+n₂*s₁)
                    ms₂₂ = (q₁*s₁+q₂*s₂)*n₂*s₂

                    W₁₁ = (q₁n₁+ms₁₁)*𝑤/2/𝐴
                    W₁₂ = (q₁n₂+q₂n₁+2*ms₁₂)*𝑤/4/𝐴
                    W₂₂ = (q₂n₂+ms₂₂)*𝑤/2/𝐴
                    for i in 1:length(𝓒)
                        ∂²𝝭∂x²[i] += 𝝭[i]*W₁₁
                        ∂²𝝭∂x∂y[i] += 𝝭[i]*W₁₂
                        ∂²𝝭∂y²[i] += 𝝭[i]*W₂₂
                    end
                end
            end
        end

        for b in Γᶿ
            if ap∩b ≠ nothing
                for ξ in b.𝓖
                    𝑤 = ξ.𝑤
                    n₁ = ξ.n₁
                    n₂ = ξ.n₂
                    s₁ = -n₂
                    s₂ = n₁
                    ∂𝝭∂x = ξ[:∂𝝭∂x]
                    ∂𝝭∂y = ξ[:∂𝝭∂y]
                    𝒒 = get𝒑₂(b,ξ)

                    𝒒̂ᵀ𝗚⁻¹𝒒 =  𝒒̂ᵀ𝗚⁻¹*𝒒

                    mn₁₁n₁ = 𝒒̂ᵀ𝗚⁻¹𝒒*n₁*n₁*n₁
                    mn₁₁n₂ = 𝒒̂ᵀ𝗚⁻¹𝒒*n₁*n₁*n₂
                    mn₁₂n₁ = 𝒒̂ᵀ𝗚⁻¹𝒒*n₁*n₂*n₁
                    mn₁₂n₂ = 𝒒̂ᵀ𝗚⁻¹𝒒*n₁*n₂*n₂
                    mn₂₂n₁ = 𝒒̂ᵀ𝗚⁻¹𝒒*n₂*n₂*n₁
                    mn₂₂n₂ = 𝒒̂ᵀ𝗚⁻¹𝒒*n₂*n₂*n₂

                    W₁₁₁ = mn₁₁n₁*𝑤
                    W₁₁₂ = mn₁₁n₂*𝑤
                    W₁₂₁ = mn₁₂n₁*𝑤
                    W₁₂₂ = mn₁₂n₂*𝑤
                    W₂₂₁ = mn₂₂n₁*𝑤
                    W₂₂₂ = mn₂₂n₂*𝑤

                    for i in 1:length(𝓒)
                        ∂²𝝭∂x²[i] += ∂𝝭∂x[i]*W₁₁₁ + ∂𝝭∂y[i]*W₁₁₂
                        ∂²𝝭∂x∂y[i] += ∂𝝭∂x[i]*W₁₂₁ + ∂𝝭∂y[i]*W₁₂₂
                        ∂²𝝭∂y²[i] += ∂𝝭∂x[i]*W₂₂₁ + ∂𝝭∂y[i]*W₂₂₂
                    end
                end
            end
        end

        for c in Γᴾ
            if ap∩c ≠ nothing
                ξ = c.𝓖[1]
                Δn₁s₁ = ξ.Δn₁s₁
                Δn₁s₂n₂s₁ = ξ.Δn₁s₂n₂s₁
                Δn₂s₂ = ξ.Δn₂s₂

                𝝭 = ξ[:𝝭]
                𝒒 = get𝒑₂(c,ξ)

                𝒒̂ᵀ𝗚⁻¹𝒒 =  𝒒̂ᵀ𝗚⁻¹*𝒒

                Δms₁₁ =  𝒒̂ᵀ𝗚⁻¹𝒒*Δn₁s₁
                Δms₂₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*Δn₂s₂
                Δms₁₂ =  𝒒̂ᵀ𝗚⁻¹𝒒*Δn₁s₂n₂s₁/2

                for i in 1:length(𝓒)
                    ∂²𝝭∂x²[i] += 𝝭[i]*Δms₁₁
                    ∂²𝝭∂x∂y[i] += 𝝭[i]*Δms₁₂
                    ∂²𝝭∂y²[i] += 𝝭[i]*Δms₂₂
                end
            end
        end
    end
end

for set𝝭 in (:set𝝭!,:set∇𝝭!,:set∇₂𝝭!,:set∇²𝝭!,:set∇³𝝭!,:set∇̂³𝝭!,:set∇²₂𝝭!)
    @eval begin
        function $set𝝭(aps::Vector{T}) where T<:ReproducingKernel
            for ap in aps
                𝓖 = ap.𝓖
                for 𝒙 in 𝓖
                    $set𝝭(ap,𝒙)
                end
            end
        end
    end
end

for set𝝭 in (:set∇̃𝝭!,:set∇̃²𝝭!,:set∇∇̃²𝝭!)
    @eval begin
        function $set𝝭(gps::Vector{T},aps::Vector{S}) where {T<:ReproducingKernel,S<:ReproducingKernel}
            if length(gps) ≠ length(aps)
                error("Miss match element numbers")
            else
                for i in 1:length(gps)
                    $set𝝭(gps[i],aps[i])
                end
            end
        end
    end
end

function set∇̄𝝭!(aps::Vector{T}) where T<:ReproducingKernel
    for ap in aps
        set∇̄𝝭!(ap)
    end
end

function set∇̄²𝝭!(aps::Vector{T};Γᵍ::Vector{T}=T[],Γᶿ::Vector{T}=T[],Γᴾ::Vector{T}=T[]) where T<:ReproducingKernel
    for ap in aps
        set∇̄²𝝭!(ap,Γᵍ=Γᵍ,Γᶿ=Γᶿ,Γᴾ=Γᴾ)
    end
end

function set∇∇̄²𝝭!(aps::Vector{T};Γᵍ::Vector{T}=T[],Γᶿ::Vector{T}=T[],Γᴾ::Vector{T}=T[]) where T<:ReproducingKernel
    for i in 1:length(aps)
        isempty(Γᵍ) ? a = nothing : a = Γᵍ[i]
        isempty(Γᶿ) ? b = nothing : b = Γᶿ[i]
        set∇∇̄²𝝭!(aps[i],Γᵍ=a,Γᶿ=b,Γᴾ=Γᴾ)
    end
end
