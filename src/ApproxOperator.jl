module ApproxOperator

import Base: +, -, *, getindex, getproperty, setproperty!, length, push!
import InteractiveUtils: subtypes

abstract type AbstractNode end
abstract type Approximator end
abstract type AbstractPoi <: Approximator end
abstract type AbstractSeg <: Approximator end
abstract type AbstractTri <: Approximator end
abstract type AbstractQuad <: Approximator end
abstract type AbstractTet <: Approximator end

include("node.jl")
export Node
include("approximation.jl")
export Poi1,Seg2
include("operation.jl")
include("import.jl")
# export import_msh
# include("export.jl")
# # export export_shape_functions, VTKExport
# include("efficiency.jl")
# export efficiency
# export efficiency, f
export Operator

end
