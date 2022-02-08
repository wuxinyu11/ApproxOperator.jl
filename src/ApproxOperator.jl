module ApproxOperator

import Base: +, -, *, getindex, setindex!, getproperty, setproperty!, length, push!, fill!, similar
import InteractiveUtils: subtypes

abstract type AbstractNode end
abstract type AbstractElement{T} end
abstract type SpatialPartition end

include("node.jl")
include("approximation.jl")
include("meshfree.jl")
include("integration.jl")
include("operation.jl")
include("import.jl")
export Node, Element, SNode, ReproducingKernel
export importmsh
export RegularGrid
export Operator, prescribe!, similar, glue
export set𝓖!
export set𝝭!, set∇𝝭!, set∇̃𝝭!, set𝒏!

end
