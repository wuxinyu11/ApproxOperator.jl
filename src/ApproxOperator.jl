module ApproxOperator

import Base: +, -, *, getindex, setindex!, getproperty, setproperty!, length, push!, fill!, issubset, intersect
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
export Node, Element, SNode, ReproducingKernel, getnₚ
export importmsh
export RegularGrid
export Operator, prescribe!, issubset, intersect
export set𝓖!
export set𝝭!, set∇𝝭!, set∇²𝝭!, set∇³𝝭!, set∇̃𝝭!, set∇̃²𝝭!, set∇̄𝝭!, set𝒏!, setg̃!, set∇𝑢!, get∇𝑢, get𝝐

#debug
include("littletools.jl")
export get𝐴,cal𝗠!,cal𝗚!,get𝒙,get∇𝝭,get𝝭,checkIC, checkCC, checkConsistency

end
