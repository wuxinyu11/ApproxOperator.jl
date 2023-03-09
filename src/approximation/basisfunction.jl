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
@inline get∇𝒑₁(ap::ReproducingKernel{:Linear1D,𝑠,𝜙,T},ξ::Any) where {𝑠,𝜙,T} = get𝒑₁(ap,ξ), get∂𝒑₁∂ξ(ap,ξ)
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
@inline get𝒑₁(ap::ReproducingKernel{:Quadratic1D},ξ::Node) = get𝒑₁(ap,ξ.ξ)
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
@inline get𝒑₁(ap::ReproducingKernel{:Cubic1D},ξ::Node) = get𝒑₁(ap,ξ.ξ)
@inline get𝒑₁(::ReproducingKernel{:Cubic1D},ξ::Float64) = (1.0,0.5*(1.0-ξ),0.25*(1.0-ξ)^2)
@inline get∂𝒑₁∂ξ(ap::ReproducingKernel{:Cubic1D},ξ::Node) = get∂𝒑₁∂ξ(ap,ξ.ξ)
@inline get∂𝒑₁∂ξ(::ReproducingKernel{:Cubic1D},ξ::Float64) = (0.,1.0,(1.0-ξ))

# ------------ Linear2D ---------------
@inline get𝑛𝒑(::ReproducingKernel{:Linear2D}) = 3
@inline get𝒑(::ReproducingKernel{:Linear2D},x::NTuple{3,Float64}) = (1.,x[1],x[2])
@inline get∂𝒑∂x(::ReproducingKernel{:Linear2D}, ::Any) = (0.,1.,0.)
@inline get∂𝒑∂y(::ReproducingKernel{:Linear2D}, ::Any) = (0.,0.,1.)
@inline get∂𝒑∂z(::ReproducingKernel{:Linear2D}, ::Any) = (0.,0.,0.)

@inline get𝑛𝒑₁(::ReproducingKernel{:Linear2D}) = 1
@inline get𝒑₁(ap::ReproducingKernel{:Linear2D},ξ::Node) = get𝒑₁(ap,ξ.ξ,ξ.η)
@inline get𝒑₁(::ReproducingKernel{:Linear2D},::Any,::Any) = (1.,)
@inline get∂𝒑₁∂ξ(ap::ReproducingKernel{:Linear2D},ξ::Node) = get∂𝒑₁∂ξ(ap,ξ.ξ,ξ.η)
@inline get∂𝒑₁∂ξ(::ReproducingKernel{:Linear2D},::Any,::Any) = (0.,)
@inline get∂𝒑₁∂η(ap::ReproducingKernel{:Linear2D},ξ::Node) = get∂𝒑₁∂η(ap,ξ.ξ,ξ.η)
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
@inline get𝒑₁(ap::ReproducingKernel{:Quadratic2D},ξ::Node) = get𝒑₁(ap,ξ.ξ,ξ.η)
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
@inline get𝒑₁(ap::ReproducingKernel{:Cubic2D},ξ::Node) = get𝒑₁(ap,ξ.ξ,ξ.η)
@inline get𝒑₁(::ReproducingKernel{:Cubic2D},ξ::Float64,η::Float64) = (1.,ξ,η,ξ^2,ξ*η,η^2)
@inline get∂𝒑₁∂ξ(ap::ReproducingKernel{:Cubic2D},ξ::Node) = get∂𝒑₁∂ξ(ap,ξ.ξ,ξ.η)
@inline get∂𝒑₁∂ξ(::ReproducingKernel{:Cubic2D},ξ::Float64,η::Float64) = (0.,1.,0.,2.0*ξ,η,0.)
@inline get∂𝒑₁∂η(ap::ReproducingKernel{:Cubic2D},ξ::Node) = get∂𝒑₁∂η(ap,ξ.ξ,ξ.η)
@inline get∂𝒑₁∂η(::ReproducingKernel{:Cubic2D},ξ::Float64,η::Float64) = (0.,0.,1.,0.,ξ,2.0*η)

@inline get𝑛𝒑₂(::ReproducingKernel{:Cubic2D}) = 3
@inline get𝒑₂(ap::ReproducingKernel{:Cubic2D},ξ::Node) = get𝒑₂(ap,ξ.ξ,ξ.η)
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
@inline get𝒑₁(ap::ReproducingKernel{:Quartic2D},ξ::Node) = get𝒑₁(ap,ξ.ξ,ξ.η)
@inline get𝒑₁(::ReproducingKernel{:Quartic2D},ξ::Float64,η::Float64) = (1.,ξ,η,ξ^2,ξ*η,η^2,ξ^3,ξ^2*η,ξ*η^2,η^3)
@inline get∂𝒑₁∂ξ(ap::ReproducingKernel{:Quartic2D},ξ::Node) = get∂𝒑₁∂ξ(ap,ξ.ξ,ξ.η)
@inline get∂𝒑₁∂ξ(::ReproducingKernel{:Quartic2D},ξ::Float64,η::Float64) = (0.,1.,0.,2.0*ξ,η,0.,3.0*ξ^2,2.0*ξ*η,η^2,0.)
@inline get∂𝒑₁∂η(ap::ReproducingKernel{:Quartic2D},ξ::Node) = get∂𝒑₁∂η(ap,ξ.ξ,ξ.η)
@inline get∂𝒑₁∂η(::ReproducingKernel{:Quartic2D},ξ::Float64,η::Float64) = (0.,0.,1.,0.,ξ,2.0*η,0.,ξ^2,2.0*ξ*η,3.0*η^2)

@inline get𝑛𝒑₂(::ReproducingKernel{:Quartic2D}) = 6
@inline get𝒑₂(ap::ReproducingKernel{:Quartic2D},ξ::Node) = get𝒑₂(ap,ξ.ξ,ξ.η)
@inline get𝒑₂(ap::ReproducingKernel{:Quartic2D},ξ::NTuple{3,Float64}) = get𝒑₂(ap,ξ[1],ξ[2])
@inline get𝒑₂(::ReproducingKernel{:Quartic2D},ξ::Float64,η::Float64) = (1.,ξ,η,ξ^2,ξ*η,η^2)
@inline get∂𝒑₂∂ξ(ap::ReproducingKernel{:Quartic2D},ξ::Node) = get∂𝒑₂∂ξ(ap,ξ.ξ,ξ.η)
@inline get∂𝒑₂∂ξ(ap::ReproducingKernel{:Quartic2D},ξ::Float64,η::Float64) = (0.,1.,0.,2.0*ξ,η,0.)
@inline get∂𝒑₂∂η(ap::ReproducingKernel{:Quartic2D},ξ::Node) = get∂𝒑₂∂η(ap,ξ.ξ,ξ.η)
@inline get∂𝒑₂∂η(ap::ReproducingKernel{:Quartic2D},ξ::Float64,η::Float64) = (0.,0.,1.,0.,ξ,2.0*η)
@inline get∂²𝒑₂∂ξ²(ap::ReproducingKernel{:Quartic2D},ξ::Any) = (0.,0.,0.,2.,0.,0.)
@inline get∂²𝒑₂∂ξ∂η(ap::ReproducingKernel{:Quartic2D},ξ::Any) = (0.,0.,0.,0.,1.,0.)
@inline get∂²𝒑₂∂η²(ap::ReproducingKernel{:Quartic2D},ξ::Any) = (0.,0.,0.,0.,0.,2.)

