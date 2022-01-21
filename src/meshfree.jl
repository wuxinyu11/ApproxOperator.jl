
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
    function (sp::t)(ap::T) where T<:Approximator
        𝓒 = ap.𝓒
        indices = Set{Int}()
        for 𝒙 in 𝓒
            union!(indices,sp(𝒙.x,𝒙.y,𝒙.z))
        end
        union!(𝓒,(Node(i,𝓒[1].data) for i in indices))
    end
    function (sp::t)(aps::Vector{T}) where T<:Approximator
        for ap in aps
            sp(ap)
        end
    end
end

## Basis Function
@inline get𝒑(ap::Approximator,x::NTuple{3,Float64}) = get𝒑(ap.type[1],x)
@inline get∂𝒑∂x(ap::Approximator,x::NTuple{3,Float64}) = get∂𝒑∂x(ap.type[1],x)
@inline get∂𝒑∂y(ap::Approximator,x::NTuple{3,Float64}) = get∂𝒑∂y(ap.type[1],x)
@inline get∂𝒑∂z(ap::Approximator,x::NTuple{3,Float64}) = get∂𝒑∂z(ap.type[1],x)
@inline get∂²𝒑∂x²(ap::Approximator,x::NTuple{3,Float64}) = get∂²𝒑∂x²(ap.type[1],x)
@inline get∂²𝒑∂y²(ap::Approximator,x::NTuple{3,Float64}) = get∂²𝒑∂y²(ap.type[1],x)
@inline get∇𝒑(ap::Approximator,x::NTuple{3,Float64}) = get𝒑(ap.type[1],x), get∂𝒑∂x(ap.type[1],x), get∂𝒑∂y(ap.type[1],x), get∂𝒑∂z(ap.type[1],x)

# ------------ Linear1D ---------------
@inline get𝒑(::Val{:Linear1D},x::NTuple{3,Float64}) = (1.,x[1])
@inline get∂𝒑∂x(::Val{:Linear1D}, ::NTuple{3,Float64}) = (0.,1.)
@inline get∂𝒑∂y(::Val{:Linear1D}, ::NTuple{3,Float64}) = (0.,0.)
@inline get∂𝒑∂z(::Val{:Linear1D}, ::NTuple{3,Float64}) = (0.,0.)

# ------------ Quadaratic1D ---------------
@inline get𝒑(::Val{:Quadratic1D},x::NTuple{3,Float64}) = (1.,x[1],x[1]^2)
@inline get∂𝒑∂x(::Val{:Quadratic1D},x::NTuple{3,Float64}) = (0.,1.,2*x[1])
@inline get∂𝒑∂y(::Val{:Quadratic1D}, ::NTuple{3,Float64}) = (0.,0.,0.)
@inline get∂𝒑∂z(::Val{:Quadratic1D}, ::NTuple{3,Float64}) = (0.,0.,0.)
@inline get∂²𝒑∂x²(::Val{:Quadratic1D}, ::NTuple{3,Float64}) =(0.,0.,2.)

# ------------ Cubic1D ---------------
@inline get𝒑(::Val{:Cubic1D},x::NTuple{3,Float64}) = (1.,x[1],x[1]^2,x[1]^3)
@inline get∂𝒑∂x(::Val{:Cubic1D},x::NTuple{3,Float64}) = (0.,1.,2*x[1],3*x[1]^2)
@inline get∂𝒑∂y(::Val{:Cubic1D}, ::NTuple{3,Float64}) = (0.,0.,0.,0.)
@inline get∂𝒑∂z(::Val{:Cubic1D}, ::NTuple{3,Float64}) = (0.,0.,0.,0.)
@inline get∂²𝒑∂x²(::Val{:Cubic1D},x::NTuple{3,Float64}) = (0.,0.,2.,6*x[1])

# ------------ Linear2D ---------------
@inline get𝒑(::Val{:Linear2D},x::NTuple{3,Float64}) = (1.,x[1],x[2])
@inline get∂𝒑∂x(::Val{:Linear2D}, ::NTuple{3,Float64}) = (0.,1.,0.)
@inline get∂𝒑∂y(::Val{:Linear2D}, ::NTuple{3,Float64}) = (0.,0.,1.)
@inline get∂𝒑∂z(::Val{:Linear2D}, ::NTuple{3,Float64}) = (0.,0.,0.)

# ------------ Quadratic2D ---------------
@inline get𝒑(::Val{:Quadratic2D},x::NTuple{3,Float64}) = (1.,x[1],x[2],x[1]^2,x[1]*x[2],x[2]^2)
@inline get∂𝒑∂x(::Val{:Quadratic2D},x::NTuple{3,Float64}) = (0.,1.,0.,2*x[1],x[2],0.)
@inline get∂𝒑∂y(::Val{:Quadratic2D},x::NTuple{3,Float64}) = (0.,0.,1.,0.,x[1],2*x[2])
@inline get∂𝒑∂z(::Val{:Quadratic2D}, ::NTuple{3,Float64}) = (0.,0.,0.,0.,0.,0.)

# ------------ Cubic2D ---------------
@inline get𝒑(::Val{:Cubic2D},x::NTuple{3,Float64}) =
(
    1., x[1], x[2], x[1]^2, x[1]*x[2], x[2]^2, x[1]^3, x[1]^2*x[2], x[1]*x[2]^2, x[2]^3
)
@inline get∂𝒑∂x(::Val{:Cubic2D},x::NTuple{3,Float64}) =
(
    0., 1., 0., 2*x[1], x[2], 0., 3*x[1]^2, 2*x[1]*x[2], x[2]^2, 0.
)
@inline get∂𝒑∂y(::Val{:Cubic2D},x::NTuple{3,Float64}) =
(
    0., 0., 1., 0., x[1], 2*x[2], 0., x[1]^2, 2*x[1]*x[2], 3*x[2]^2
)
@inline get∂𝒑∂z(::Val{:Cubic2D},::NTuple{3,Float64}) =
(
    0., 0., 0., 0., 0., 0., 0., 0., 0., 0.
)

## Kernel Function
@inline get𝜙(ap::Approximator,xᵢ::Node,x::NTuple{3,Float64}) = get𝜙(ap.type[2],ap.type[3],xᵢ,x)
@inline get∂𝜙∂x(ap::Approximator,xᵢ::Node,x::NTuple{3,Float64}) = get∂𝜙∂x(ap.type[2],ap.type[3],xᵢ,x)
@inline get∇𝜙(ap::Approximator,xᵢ::Node,x::NTuple{3,Float64}) = get∇𝜙(ap.type[2],ap.type[3],xᵢ,x)

function get𝜙(::Val{:□},type::Val,x::Node,Δx::NTuple{3,Float64})
    rx = abs(Δx[1])/x.s₁
    ry = abs(Δx[2])/x.s₁
    rz = abs(Δx[3])/x.s₁
    wx = get𝜙ᵣ(type,rx)
    wy = get𝜙ᵣ(type,ry)
    wz = get𝜙ᵣ(type,rz)
    return wx*wy*wz
end

function get∂𝜙∂x(::Val{:□},type::Val,x::Node,Δx::NTuple{3,Float64})
    rx = abs(Δx[1])/x.s₁
    ∂rx = sign(Δx[1])/x.s₁
    wx = get𝜙ᵣ(type,rx)
    ∂wx = get∂𝜙∂r(type,rx)*∂rx
    return wx, ∂wx
end

function get∇𝜙(::Val{:□},type::Val,x::Node,Δx::NTuple{3,Float64})
    rx = abs(Δx[1])/x.s₁
    ry = abs(Δx[2])/x.s₂
    rz = abs(Δx[3])/x.s₃
    ∂rx = sign(Δx[1])/x.s₁
    ∂ry = sign(Δx[2])/x.s₂
    ∂rz = sign(Δx[3])/x.s₃
    wx = get𝜙ᵣ(type,rx)
    wy = get𝜙ᵣ(type,ry)
    wz = get𝜙ᵣ(type,rz)
    ∂wx = get∂𝜙∂r(type,rx)*∂rx
    ∂wy = get∂𝜙∂r(type,ry)*∂ry
    ∂wz = get∂𝜙∂r(type,rz)*∂rz
    return wx*wy*wz, ∂wx*wy*wz, wx*∂wy*wz, wx*wy*∂wz
end

## --------------- Kernel ---------------
function get𝜙ᵣ(::Val{:CubicSpline},r::Float64)
    if r > 1.0
        return 0.0
    elseif r <= 0.5
        return 2/3 - 4*r^2 +  4*r^3
    else
        return 4/3 - 4*r + 4*r^2 - 4*r^3/3
    end
end

function get∂𝜙∂r(::Val{:CubicSpline},r::Float64)
    if r > 1.0
        return 0.0
    elseif r <= 0.5
        return - 8*r + 12*r^2
    else
        return - 4   + 8*r - 4*r^2
    end
end

function get∂²𝜙∂r²(::Val{:CubicSpline},r::Float64)
    if r > 1.0
        return 0.0
    elseif r <= 0.5
        return - 8 + 24*r
    else
        return   8 - 8*r
    end
end

## calulate shape functions
function cal𝗠!(ap::T,x::NTuple{3,Float64}) where T<:ReproducingKernel
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

function cal𝗚!(dp::A) where A<:ReproducingKernel
    𝓖 = dp.𝓖
    𝗚 = ap.𝗚[:∂1]
    n = length(get𝒑(dp.type),(0.0,0.0,0.0))
    fill!(𝗚,0.0)
    for ξ in 𝓖
        w = ξ.w
        𝒑 = get𝒑(dp,ξ)
        for I in 1:n
            for J in 1:n
                𝗚[I,J] += w*𝒑[I]*𝒑[J]
            end
        end
    end
    cholesky!(𝗚)
    U⁻¹ = inverse!(𝗚)
    𝗚⁻¹ = UUᵀ!(U⁻¹)
    return 𝗚⁻¹
end

## Shape functions
function (op::Operator{:𝝭})(ap::ReproducingKernel{SNode})
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

function (op::Operator{:∇𝝭})(ap::ReproducingKernel{SNode})
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

function Operator(t::Val{:𝝭ʳ})
    id = Dict{NTuple{3,Float64},Tuple{Int,Int}}()
    ids = Int[]
    index = Int[0]
    return Operator(t,Dict(:n=>0,:id=>id,:ids=>ids,:index=>index))
end

function (op::Operator{:𝝭ʳ})(ap::ReproducingKernel{SNode})
    𝓒 = ap.𝓒
    𝓖 = ap.𝓖
    for ξ in 𝓖
        x = getx(ap,ξ)
        if haskey(op.id,x)
            i,n = op.id[x]
            index = op.index[n]+1:op.index[n+1]
            ids = @views op.ids[index]
            for j in 1:length(𝓒)
                id = findfirst(x->x==𝓒[j].id,ids)
                if id ≠ nothing
                    ξ.𝝭[:∂1][ξ.index[ξ.id]+j] = ξ.𝝭[:∂1][ξ.index[i]+id]
                end
            end
        else
            i = ξ.id
            I = ξ.index[i]
            ξ̂ = Node(ξ)
            𝝭 = get𝝭(ap,ξ̂)
            for j in 1:length(𝓒)
                ξ.𝝭[:∂1][I+j] = 𝝭[j]
            end
            op.n += 1
            push!(op.id,x=>(i,op.n))
            push!(op.ids,(ξ_.id for ξ_ in 𝓒)...)
            push!(op.index,last(op.index)+length(𝓒))
        end
    end
end

## RK gradient smoothing
function Operator(t::Val{:∇̃𝝭},𝒑::Val)
    n = length(get𝒑(𝒑,(0.0,0.0,0.0)))
    𝗚 = Dict(:∂x=>SymMat(n),:∂y=>SymMat(n),:∂z=>SymMat(n))
    return Operator(t,Dict(:n=>0,:id=>id,:ids=>ids,:index=>index))
end

function (op::Operator{:∇̃𝝭})(ap::SegN{SNode})
    L = ap.L
    𝗚⁻¹ = cal𝗚!(dp)
    𝓒 = ap.𝒞
    𝓖 = SNode(dp)
    for ξ̂ in 𝓖
        ξ̂.index[ξ̂.id+1] = ξ̂.index[ξ̂.id]+length(𝓒)
        𝒑̂ = get𝒑(dp,ξ̂)
        𝒑̂ᵀ𝗚⁻¹ = 𝒑̂*𝗚⁻¹
        ∂𝝭∂x = dp.𝝭[:∂x]
        fill!(∂𝝭∂x,0.0)
        for ξ in ap.𝓖
            w = ξ.w
            wᵇ = ξ.wᵇ
            n₁ = ξ.n₁
            𝝭 = get𝝭(ap,ξ)
            𝒑, ∂𝒑∂ξ = get∇𝒑(dp,ξ)
            𝒑̂ᵀ𝗚⁻¹𝒑 = 𝒑̂ᵀ𝗚⁻¹*𝒑
            𝒑̂ᵀ𝗚⁻¹∂𝒑∂ξ = 𝒑̂ᵀ𝗚⁻¹*∂𝒑∂ξ
            for i in 1:length(𝓒)
                ∂𝝭∂x[i] += 𝝭[i]*𝒑̂ᵀ𝗚⁻¹𝒑*n₁*wᵇ + 𝝭[i]*𝒑̂ᵀ𝗚⁻¹∂𝒑∂ξ/L*w
            end
        end
        for i in 1:length(𝓒)
            ξ̂.𝝭[:∂x][ξ̂.index[ξ̂.id]+i] = ∂𝝭∂x[i]
        end
    end
    return Seg2(𝓒,𝓖)
end
