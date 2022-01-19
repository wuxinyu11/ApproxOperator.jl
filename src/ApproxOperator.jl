module ApproxOperator

import Base: +, -, *, getindex, setindex!, getproperty, setproperty!, length, push!, fill!
import InteractiveUtils: subtypes

abstract type AbstractNode end
abstract type Approximator end
abstract type AbstractPoi <: Approximator end
abstract type AbstractSeg <: Approximator end
abstract type AbstractTri <: Approximator end
abstract type AbstractQuad <: Approximator end
abstract type AbstractTet <: Approximator end
abstract type ReproducingKernel{T,𝒑,𝑠,𝜙} <: Approximator  where T<:AbstractNode end
abstract type SpatialPartition end

include("node.jl")
export Node, SNode
include("meshfree.jl")
include("approximation.jl")
export Poi1,Seg2,SegN
include("operation.jl")
include("import.jl")
# export import_msh
# include("export.jl")
# # export export_shape_functions, VTKExport
# include("efficiency.jl")
# export efficiency
# export efficiency, f
export Operator, prescribe!
export get𝒑, SymMat
end
