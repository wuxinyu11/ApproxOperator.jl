module ApproxOperator

import Base: +, -, *, getindex, length
import InteractiveUtils: subtypes

abstract type AbstractNode end
abstract type PhysicalNode<:AbstractNode end
abstract type ParametricNode<:AbstractNode end
abstract type Approximator end
abstract type BasisFunction end
abstract type KernelFunction end
abstract type ShapeFunction end
abstract type SpatialPartition end
abstract type AbstractPoi <: Approximator end
abstract type AbstractSeg <: Approximator end
abstract type AbstractTri <: Approximator end
abstract type AbstractQuad <: Approximator end
abstract type AbstractTet <: Approximator end
abstract type Operator end

include("node.jl")
export Node
include("approximation.jl")
include("operation.jl")
include("import.jl")
export import_msh
include("export.jl")
# export export_shape_functions, VTKExport
include("efficiency.jl")
export efficiency

## Meshfree
import Base: setindex!, fill!

abstract type MeshfreeSpace end
include("../meshfree/node.jl")
include("../meshfree/approximation.jl")
include("../meshfree/efficiency_meshfree.jl")

# debug
export efficiency_meshfree
export Node, Seg2, MFPoint
export RegularGrid, TensorProductKernel, MFSpace
# export get_jacobe, get_coordinates, get_shape_functions

end
