module ApproxOperator

import Base: +, -, *, /, getindex, setindex!, getproperty, setproperty!, length, push!, fill!, issubset, intersect
import InteractiveUtils: subtypes

abstract type AbstractNode end
abstract type AbstractElement{T} end
abstract type SpatialPartition end

include("node.jl")
include("element.jl")
include("meshfree.jl")
include("integration.jl")
include("operation.jl")
include("import.jl")

export getnₚ, prescribe!, Operator
export set𝓖!
export set𝝭!, set𝝭𝑛!, set𝝭̄!, set𝝭ₕ!, set∇𝝭!, set∇₁𝝭!, set∇₂𝝭!, set∇𝝭ₕ!, set∇²𝝭!, set∇²₂𝝭!, set∇³𝝭!, set∇̃𝝭!, set∇̃₁𝝭!, set∇̃²𝝭!, set∇∇̃²𝝭!, set∇̄𝝭!, set∇̄²𝝭!, set∇∇̄²𝝭!, set∇̂³𝝭!, set𝒏!, set∇𝑢!, get∇𝑢, get𝝐, set_memory_𝗠!, set_memory_𝝭!

# debug
include("littletools.jl")
export check𝝭, check∇𝝭, check∇₁𝝭, check∇₂𝝭, check∇²𝝭, check∇²₂𝝭, check∇³𝝭, check∇∇²𝝭

end