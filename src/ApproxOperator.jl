module ApproxOperator

import Base: +, -, *, getindex, setindex!, getproperty, setproperty!, length, push!, fill!, issubset, intersect
import InteractiveUtils: subtypes
# import Base: AbstractMatrix, size, strides
import LinearAlgebra.LAPACK: potrf!, potri!, sytrf!, sytri!

abstract type AbstractNode end
abstract type AbstractElement{T} end
abstract type SpatialPartition end

include("node.jl")
include("snode.jl")
include("approximation.jl")
include("approximation_mf.jl")
include("integration.jl")
include("operation.jl")
include("operation_thin_plate.jl")
include("approximation_rk.jl")
include("import.jl")

export Node, Element, SNode, ReproducingKernel, getnₚ
export importmsh
export RegularGrid
export Operator, prescribe!, issubset, intersect
export set𝓖!
export set𝝭!, set∇𝝭!, set∇²𝝭!, set∇³𝝭!, set∇̃𝝭!, set∇̃²𝝭!, set∇∇̃²𝝭!, set∇̄𝝭!, set∇̄²𝝭!, set∇∇̄²𝝭!, set𝒏!, set∇𝑢!, get∇𝑢, get𝝐, set_memory_𝝭!

#debug
include("littletools.jl")
export get𝐴,cal𝗠!,cal𝗚!,get𝒙,get∇𝝭,get𝝭,checkIC, checkCC, checkConsistency

end
