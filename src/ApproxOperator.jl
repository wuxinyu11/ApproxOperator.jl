module ApproxOperator

import Base: +, -, *, getindex, setindex!, getproperty, setproperty!, length, push!, fill!, similar
import InteractiveUtils: subtypes

abstract type AbstractNode end
abstract type Approximator end
abstract type AbstractPoi<:Approximator end
abstract type AbstractSeg<:Approximator end
abstract type AbstractTri<:Approximator end
abstract type AbstractQuad<:Approximator end
abstract type AbstractTet<:Approximator end
abstract type ReproducingKernel{T<:AbstractNode,𝒑,𝑠,𝜙} <: Approximator end
abstract type SpatialPartition end

include("operation.jl")
include("node.jl")
include("meshfree.jl")
include("approximation.jl")
include("import.jl")
export Node, SNode
export Poi1, Seg2, Tri3, Quad, PoiN, SegN
export importdata
export RegularGrid
export Operator, prescribe!, similar
export set𝓖!
export set𝝭!, set∇𝝭!, set∇̃𝝭!, set𝒏!

end
